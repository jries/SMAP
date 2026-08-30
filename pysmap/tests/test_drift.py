"""Drift correction: recover a known drift, and correct everything with it."""
import numpy as np
import pytest

from smapfit.drift import (Drift, DriftSettings, correct_drift, drift_corrected_path,
                           estimate_drift, load_drift, save_drift_corrected)
from smapfit.filter import LocFilter
from smapfit.locs import Localizations

comet = pytest.importorskip("comet")


def simulate(n_points=200, n_frames=50, per_frame=30, seed=1):
    """Localizations of a fixed random structure, moved by a known drift."""
    rng = np.random.default_rng(seed)
    structure = rng.random((n_points, 3)) * 1000.0
    t = np.arange(n_frames)
    base = (3 * np.sin(t / 20.0 + 0.37) + 1.25 * np.cos(t / 6.0 + 1.2) + t / 25) * 10.0
    drift = np.column_stack([base, 0.7 * base, -0.5 * base])

    frame = np.repeat(t, per_frame)
    rng.shuffle(frame)
    true_xyz = structure[rng.integers(0, n_points, frame.size)]
    xyz = true_xyz + drift[frame]
    locs = Localizations({"frame": frame.astype(np.int64),
                          "x_nm": xyz[:, 0].astype(np.float32),
                          "y_nm": xyz[:, 1].astype(np.float32),
                          "z_nm": xyz[:, 2].astype(np.float32),
                          "loc_precision_nm": rng.random(frame.size).astype(np.float32) * 40},
                         {"units": "nm"})
    return locs, drift, true_xyz


# the per-window fit; the package default is the spline, tested separately
SETTINGS = DriftSettings(segmentation_var=2, initial_sigma_nm=120, max_drift_nm=100,
                         target_sigma_nm=10, backend="cpu", spline=False, group=False)


def test_recovers_known_drift():
    locs, truth, _ = simulate()
    drift = estimate_drift(locs, SETTINGS)

    assert len(drift) == 50
    # both are defined up to a constant offset: compare after centring
    error = (drift.drift - drift.drift.mean(0)) - (truth - truth.mean(0))
    assert np.abs(error).mean() < 1.0


def test_corrects_localizations_the_filter_hides():
    locs, _, true_xyz = simulate()
    # estimate from the good half only, correct all of it
    keep = LocFilter(locs, loc_precision_nm=(None, 20))
    corrected, drift = correct_drift(locs, SETTINGS, select=keep)

    assert drift.n_used == len(keep.indices) < len(locs)
    assert len(corrected) == len(locs)
    # the hidden localizations must land back on the undrifted structure
    # (up to the constant offset drift is defined against)
    hidden = ~keep.mask
    for axis, column in enumerate(("x_nm", "y_nm", "z_nm")):
        residual = corrected[column][hidden] - true_xyz[hidden, axis]
        assert residual.std() < 1.0
        assert (locs[column][hidden] - true_xyz[hidden, axis]).std() > 5


def test_pixel_table_uses_the_pixel_size():
    locs, _, true_xyz = simulate()
    pixels = Localizations({"frame": locs["frame"],
                            "x_pix": locs["x_nm"] / 100.0,
                            "y_pix": locs["y_nm"] / 100.0},
                           {"units": "pixel", "pixelsize_nm": 100.0})
    drift = estimate_drift(locs, SETTINGS)
    corrected = drift.apply(pixels)

    assert np.allclose(corrected["x_pix"] * 100.0,
                       drift.apply(locs)["x_nm"], atol=1e-2)


def test_save_and_reload(tmp_path):
    locs, _, true_xyz = simulate()
    drift = Drift(np.arange(150, dtype=float).reshape(50, 3))
    path = save_drift_corrected(tmp_path / "run_driftc.hdf5", drift.apply(locs), drift)

    from smapfit.io.hdf5 import load_localizations
    assert len(load_localizations(path)) == len(locs)
    assert np.allclose(load_drift(path).drift, drift.drift)
    assert drift_corrected_path(tmp_path / "run.hdf5").name == "run_driftc.hdf5"


def test_drift_from_another_acquisition_is_refused():
    locs, _, true_xyz = simulate()
    with pytest.raises(ValueError, match="different acquisition"):
        Drift(np.zeros((10, 3))).apply(locs)


def test_grouped_estimate_still_corrects_every_localization():
    locs, truth, true_xyz = simulate()
    # two localizations of the same emitter in consecutive frames, so grouping
    # has something to collapse
    settings = DriftSettings(segmentation_var=2, initial_sigma_nm=120, max_drift_nm=100,
                             target_sigma_nm=10, backend="cpu", spline=False,
                             group=True, group_dx_nm=30.0)
    corrected, drift = correct_drift(locs, settings)

    assert drift.n_used < len(locs)          # grouping collapsed rows
    assert len(corrected) == len(locs)       # every row is still corrected
    error = (drift.drift - drift.drift.mean(0)) - (truth - truth.mean(0))
    assert np.abs(error).mean() < 5.0


def test_two_stage_matches_the_single_pass():
    locs, truth, _ = simulate()
    single = estimate_drift(locs, SETTINGS)
    two = estimate_drift(locs, DriftSettings(
        segmentation_var=2, backend="cpu", spline=False, two_stage=True,
        group_dx_nm=30.0, two_stage_radius_nm=60.0))

    for drift in (single, two):
        error = (drift.drift - drift.drift.mean(0)) - (truth - truth.mean(0))
        assert np.abs(error).mean() < 5.0


def test_spline_drift_needs_no_time_windows():
    locs, truth, _ = simulate()
    drift = estimate_drift(locs, DriftSettings(
        backend="cpu", spline=True, group=False, spline_knot_frames=5,
        max_drift_nm=100, initial_sigma_nm=120, target_sigma_nm=10))

    error = (drift.drift - drift.drift.mean(0)) - (truth - truth.mean(0))
    assert np.abs(error).mean() < 1.0
    # the curve is the fit itself, evaluated per frame -- not interpolated
    assert len(drift) == 50


def test_spline_penalty_smooths():
    locs, _, _ = simulate()
    rough, smooth = (estimate_drift(locs, DriftSettings(
        backend="cpu", spline=True, group=False, spline_knot_frames=2,
        max_drift_nm=100, initial_sigma_nm=120, target_sigma_nm=10,
        spline_penalty=p)).drift
        for p in (0.0, 1e4))
    curvature = lambda d: np.abs(np.diff(d, 2, axis=0)).mean()
    assert curvature(smooth) < curvature(rough)
