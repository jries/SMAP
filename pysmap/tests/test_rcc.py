"""Redundant cross-correlation: recover a known drift from images."""
import numpy as np
import pytest

from smapfit.rcc import RCCSettings, estimate_drift_rcc, _solve
from tests.test_drift import simulate


def test_recovers_known_drift():
    locs, truth, _ = simulate(n_points=300, n_frames=200, per_frame=200)
    drift = estimate_drift_rcc(locs, RCCSettings(
        n_timepoints=10, pixelsize_nm=10.0, max_drift_nm=300.0,
        z_pixelsize_nm=10.0, tile_nm=250.0, axial_max_drift_nm=300.0))

    assert len(drift) == 200
    error = (drift.drift - drift.drift.mean(0)) - (truth - truth.mean(0))
    assert np.abs(error[:, :2]).mean() < 15.0   # image-based: a fraction of a pixel
    assert np.abs(error[:, 2]).mean() < 15.0    # z, from sliced 1-D correlations


def test_the_redundant_solve_outvotes_a_bad_pair():
    true = np.array([0.0, 10.0, 25.0, 5.0, -12.0])
    true -= true.mean()
    shifts = true[:, None] - true[None, :]
    shifts[0, 3] = 400.0                  # one correlation peak found in the wrong place
    shifts[3, 0] = -400.0

    assert np.abs(_solve(shifts, 5) - true).max() < 1.0


def test_every_pair_is_used():
    # with only consecutive pairs the middle window would be unconstrained
    true = np.array([0.0, 3.0, -6.0, 2.0]); true -= true.mean()
    shifts = true[:, None] - true[None, :]
    assert np.allclose(_solve(shifts, 4), true, atol=1e-6)
