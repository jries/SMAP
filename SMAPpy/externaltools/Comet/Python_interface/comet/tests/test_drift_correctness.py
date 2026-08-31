"""End-to-end numerical correctness of the drift estimate.

These run on the CPU backend so they work without a GPU. The GPU backends are
covered by the equivalence tests in test_backend_equivalence.py.
"""

import numpy as np
import pytest

from comet.core.drift_optimizer import comet_run_kd
from comet.tests.conftest import drift_residual_nm, make_drifting_dataset

# Recovered drift must agree with ground truth to well under a nanometre on
# noise-free synthetic data; observed residual is ~0.06 nm.
TOLERANCE_NM = 0.5

RUN_KWARGS = dict(
    segmentation_mode=2,
    segmentation_var=2,
    initial_sigma_nm=120,
    max_drift_nm=100,
    target_sigma_nm=10,
    boxcar_width=1,
    interpolation_method="cubic",
    mode="cpu",
    display=False,
)


def test_recovers_known_drift(drifting_dataset):
    locs, gt_drift = drifting_dataset

    drift, corrected = comet_run_kd(dataset=locs.copy(), return_corrected_locs=True, **RUN_KWARGS)

    assert drift_residual_nm(drift, gt_drift) < TOLERANCE_NM
    assert corrected.shape == locs.shape


def test_output_shapes_and_frame_column(drifting_dataset):
    locs, _ = drifting_dataset
    n_frames = int(locs[:, 3].max()) + 1

    drift = comet_run_kd(dataset=locs.copy(), **RUN_KWARGS)

    # (F, 4): dx, dy, dz, frame -- one row per frame, not per segment
    assert drift.shape == (n_frames, 4)
    np.testing.assert_array_equal(drift[:, 3], np.arange(n_frames))


def test_no_nans_in_output(drifting_dataset):
    locs, _ = drifting_dataset
    drift = comet_run_kd(dataset=locs.copy(), **RUN_KWARGS)
    assert np.isfinite(drift).all()


def test_drift_free_data_yields_near_zero_drift():
    # a static sample must not invent motion
    locs, _ = make_drifting_dataset(drift_amplitude_nm=0.0)

    drift = comet_run_kd(dataset=locs.copy(), **RUN_KWARGS)

    centred = drift[:, :3] - drift[:, :3].mean(axis=0)
    assert np.abs(centred).max() < TOLERANCE_NM


def test_correction_moves_locs_towards_their_undrifted_positions(drifting_dataset):
    """Corrected localizations must sit closer to their drift-free positions.

    The global spread cannot be used here: the template molecules span a
    1000 nm cube, which swamps a ~60 nm drift. Compare against the known
    undrifted coordinates instead, up to the constant offset that drift is
    only ever determined to.
    """
    locs, gt_drift = drifting_dataset
    frame = locs[:, 3].astype(int)
    undrifted = locs[:, :3] - gt_drift[frame]

    _, corrected = comet_run_kd(dataset=locs.copy(), return_corrected_locs=True, **RUN_KWARGS)

    def offset_free_rms(a, b):
        residual = a - b
        return float(np.sqrt(((residual - residual.mean(axis=0)) ** 2).sum(axis=1).mean()))

    before = offset_free_rms(locs[:, :3], undrifted)
    after = offset_free_rms(corrected[:, :3], undrifted)

    assert after < 0.1 * before, "correction left {:.2f} nm of {:.2f} nm".format(after, before)


def test_2d_dataset_with_zero_z_column():
    """2D data is passed in as a (N, 4) array with z all zeros.

    Uses a smaller dataset than the other cases: collapsing z makes every
    point coplanar, which multiplies the number of neighbour pairs within the
    search radius and makes the CPU backend far slower.
    """
    locs, gt_drift = make_drifting_dataset(n_points=100, n_frames=30, locs_per_frame=20)
    locs[:, 2] = 0.0
    gt_2d = gt_drift.copy()
    gt_2d[:, 2] = 0.0

    drift = comet_run_kd(dataset=locs.copy(), **RUN_KWARGS)

    # x and y are still recovered, and z must stay flat
    assert drift_residual_nm(drift[:, :2], gt_2d[:, :2]) < TOLERANCE_NM
    assert np.abs(drift[:, 2] - drift[:, 2].mean()).max() < TOLERANCE_NM


@pytest.mark.parametrize("interpolation", ["cubic", "catmull-rom"])
def test_interpolation_methods_agree(drifting_dataset, interpolation):
    locs, gt_drift = drifting_dataset
    kwargs = dict(RUN_KWARGS, interpolation_method=interpolation)

    drift = comet_run_kd(dataset=locs.copy(), **kwargs)

    # catmull-rom leaves the first and last spans undefined (zero), so compare
    # only the interior where both methods are defined
    interior = slice(len(drift) // 5, -len(drift) // 5)
    assert drift_residual_nm(drift[interior], gt_drift[interior]) < TOLERANCE_NM


@pytest.mark.parametrize("segmentation_mode,segmentation_var", [(0, 25), (1, 60), (2, 2)])
def test_all_segmentation_modes_recover_drift(drifting_dataset, segmentation_mode, segmentation_var):
    locs, gt_drift = drifting_dataset
    kwargs = dict(RUN_KWARGS, segmentation_mode=segmentation_mode,
                  segmentation_var=segmentation_var)

    drift = comet_run_kd(dataset=locs.copy(), **kwargs)

    assert drift_residual_nm(drift, gt_drift) < TOLERANCE_NM


def test_localization_noise_is_tolerated():
    locs, gt_drift = make_drifting_dataset(jitter_nm=5.0)

    drift = comet_run_kd(dataset=locs.copy(), **RUN_KWARGS)

    # looser bound: 5 nm per-localization noise propagates into the estimate
    assert drift_residual_nm(drift, gt_drift) < 3.0


class TestSaveCorrectedMoleculeSet:
    """save_corrected_locs=True must write the *corrected* locs with real frames.

    It previously wrote `sorted_dataset`, the segmented and possibly downsampled
    copy taken before correction, whose last column holds segment ids.
    """

    @staticmethod
    def run_and_load(tmp_path, locs, **extra):
        from comet.core.io_utils import load_normal_molecule_set

        out = tmp_path / "corrected.h5"
        drift, corrected = comet_run_kd(
            dataset=locs.copy(), return_corrected_locs=True,
            save_corrected_locs=True, save_filepath=str(out), **dict(RUN_KWARGS, **extra))
        return drift, corrected, load_normal_molecule_set(str(out))

    def test_saved_locs_are_drift_corrected(self, tmp_path, drifting_dataset):
        locs, gt_drift = drifting_dataset
        frame = locs[:, 3].astype(int)
        undrifted = locs[:, :3] - gt_drift[frame]

        def offset_free_rms(a, b):
            r = a - b
            return float(np.sqrt(((r - r.mean(axis=0)) ** 2).sum(axis=1).mean()))

        _, corrected, saved = self.run_and_load(tmp_path, locs)

        before = offset_free_rms(locs[:, :3], undrifted)
        after = offset_free_rms(saved[:, :3], undrifted)
        assert after < 0.1 * before, (
            "saved molecule set is not drift corrected: {:.2f} nm of {:.2f} nm "
            "remaining".format(after, before))

    def test_saved_frame_numbers_are_frames_not_segment_ids(self, tmp_path, drifting_dataset):
        locs, _ = drifting_dataset
        _, _, saved = self.run_and_load(tmp_path, locs)

        np.testing.assert_array_equal(
            np.sort(saved[:, 3].astype(int)), np.sort(locs[:, 3].astype(int)))

    def test_all_localizations_are_saved(self, tmp_path, drifting_dataset):
        locs, _ = drifting_dataset
        _, _, saved = self.run_and_load(tmp_path, locs)
        assert len(saved) == len(locs)

    def test_pixel_size_is_recorded(self, tmp_path, drifting_dataset):
        """Molecule sets store pixels, so the saved pixel size must be the one asked for."""
        import h5py

        locs, _ = drifting_dataset
        out = tmp_path / "px.h5"
        comet_run_kd(dataset=locs.copy(), save_corrected_locs=True, save_filepath=str(out),
                     pixelsize_nm=117.0, pixelsize_z_nm=250.0, **RUN_KWARGS)

        with h5py.File(str(out), "r") as f:
            assert f["molecule_set_data"]["xy_pixel_size_um"][()] == pytest.approx(0.117)
            assert f["molecule_set_data"]["z_pixel_size_um"][()] == pytest.approx(0.250)

    def test_round_trip_is_independent_of_pixel_size(self, tmp_path, drifting_dataset):
        """Positions are stored in pixels, so any pixel size must round-trip."""
        locs, _ = drifting_dataset

        _, corrected_a, saved_a = self.run_and_load(tmp_path, locs, pixelsize_nm=100.0)
        np.testing.assert_allclose(saved_a[:, :3], corrected_a[:, :3], rtol=1e-4, atol=1e-2)


def test_too_coarse_segmentation_reports_a_clear_error(drifting_dataset):
    """More frames per window than the movie has frames must not fail obscurely."""
    locs, _ = drifting_dataset          # 50 frames
    kwargs = dict(RUN_KWARGS, segmentation_mode=2, segmentation_var=10 ** 4)

    with pytest.raises(ValueError, match="segmentation parameter is likely too coarse"):
        comet_run_kd(dataset=locs.copy(), **kwargs)


def test_default_backend_works_without_a_gpu(drifting_dataset):
    """Omitting `mode` must work on any machine.

    The documented Python API example passes no `mode`, so a cuda default would
    make the README's first code block crash for most users.
    """
    locs, gt_drift = drifting_dataset
    kwargs = {k: v for k, v in RUN_KWARGS.items() if k != "mode"}

    drift = comet_run_kd(dataset=locs.copy(), **kwargs)

    assert np.isfinite(drift).all()
    assert drift_residual_nm(drift, gt_drift) < 1.0


def test_library_call_is_silent_by_default(drifting_dataset, capsys):
    """A library must not print progress unless display=True."""
    locs, _ = drifting_dataset

    comet_run_kd(dataset=locs.copy(), **RUN_KWARGS)

    assert capsys.readouterr().out == ""
