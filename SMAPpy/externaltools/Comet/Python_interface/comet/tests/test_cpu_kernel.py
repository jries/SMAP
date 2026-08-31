"""Tests for the compiled CPU cost kernel.

The plain-Python reference implementation in cpu_wrapper is the specification;
the compiled kernels must agree with it, and the analytic gradient must agree
with finite differences.
"""

import numpy as np
import pytest

from comet.core.cpu_wrapper import (
    PARALLEL_PAIR_THRESHOLD,
    _cost_and_gradient_njit,
    _cost_and_gradient_njit_parallel,
    _cost_and_gradient_reference,
    cpu_wrapper_chunked,
)
from comet.core.pair_indices import pair_indices_kdtree
from comet.core.segmenter import segmentation_wrapper
from comet.tests.conftest import make_drifting_dataset

SIGMA = 100.0
SIGMA_FACTOR = 1.0


def build_problem(n_points=300, n_frames=60, locs_per_frame=30, radius=100.0, seed=0):
    """A realistic (coords, times, idx_i, idx_j, mu) tuple straight from the pipeline."""
    locs, _ = make_drifting_dataset(n_points=n_points, n_frames=n_frames,
                                    locs_per_frame=locs_per_frame)
    result = segmentation_wrapper(locs[:, -1], 2, 2, None, return_param_dict=True)
    sorted_ds = locs.copy()
    sorted_ds[:, -1] = result.loc_segments
    sorted_ds = sorted_ds[result.loc_valid]

    idx_i, idx_j, ok = pair_indices_kdtree(sorted_ds[:, :3], radius)
    assert ok

    coords = np.ascontiguousarray(sorted_ds[:, :3], dtype=np.float64)
    times = np.ascontiguousarray(sorted_ds[:, 3], dtype=np.int64)
    idx_i = np.ascontiguousarray(idx_i, dtype=np.int64)
    idx_j = np.ascontiguousarray(idx_j, dtype=np.int64)
    # evaluate away from mu = 0, so a gradient bug cannot hide at a symmetric point
    mu = np.random.default_rng(seed).normal(0.0, 5.0, (result.n_segments, 3))
    return coords, times, idx_i, idx_j, mu


@pytest.fixture(scope="module")
def problem():
    return build_problem()


@pytest.fixture(scope="module")
def reference(problem):
    coords, times, idx_i, idx_j, mu = problem
    return _cost_and_gradient_reference(coords, times, idx_i, idx_j, mu, SIGMA, SIGMA_FACTOR)


class TestAgreesWithReference:
    def test_serial_cost(self, problem, reference):
        coords, times, idx_i, idx_j, mu = problem
        total, _ = _cost_and_gradient_njit(coords, times, idx_i, idx_j, mu, SIGMA, SIGMA_FACTOR)
        assert total == pytest.approx(reference[0], rel=1e-12)

    def test_serial_gradient(self, problem, reference):
        coords, times, idx_i, idx_j, mu = problem
        _, deri = _cost_and_gradient_njit(coords, times, idx_i, idx_j, mu, SIGMA, SIGMA_FACTOR)
        np.testing.assert_allclose(deri, reference[1], rtol=1e-10, atol=1e-12)

    def test_parallel_cost(self, problem, reference):
        coords, times, idx_i, idx_j, mu = problem
        total, _ = _cost_and_gradient_njit_parallel(coords, times, idx_i, idx_j, mu,
                                                    SIGMA, SIGMA_FACTOR)
        # different summation order, so only float-associativity slack is allowed
        assert total == pytest.approx(reference[0], rel=1e-10)

    def test_parallel_gradient(self, problem, reference):
        coords, times, idx_i, idx_j, mu = problem
        _, deri = _cost_and_gradient_njit_parallel(coords, times, idx_i, idx_j, mu,
                                                   SIGMA, SIGMA_FACTOR)
        np.testing.assert_allclose(deri, reference[1], rtol=1e-9, atol=1e-11)

    @pytest.mark.parametrize("sigma_factor", [0.1, 0.5, 1.0, 3.0])
    def test_across_the_sigma_schedule(self, problem, sigma_factor):
        """The optimizer sweeps sigma_factor down; agreement must hold throughout."""
        coords, times, idx_i, idx_j, mu = problem
        ref = _cost_and_gradient_reference(coords, times, idx_i, idx_j, mu, SIGMA, sigma_factor)
        got = _cost_and_gradient_njit(coords, times, idx_i, idx_j, mu, SIGMA, sigma_factor)

        assert got[0] == pytest.approx(ref[0], rel=1e-10)
        np.testing.assert_allclose(got[1], ref[1], rtol=1e-9, atol=1e-12)


class TestGradientIsCorrect:
    """Check the analytic gradient against finite differences.

    This is independent of the reference implementation: it would catch a
    gradient that is wrong in both.
    """

    def test_matches_finite_differences(self, problem):
        coords, times, idx_i, idx_j, mu = problem

        def cost(m):
            return _cost_and_gradient_njit(coords, times, idx_i, idx_j, m, SIGMA, SIGMA_FACTOR)[0]

        _, analytic = _cost_and_gradient_njit(coords, times, idx_i, idx_j, mu, SIGMA, SIGMA_FACTOR)

        rng = np.random.default_rng(1)
        step = 1e-6
        for _ in range(12):
            seg = int(rng.integers(0, mu.shape[0]))
            dim = int(rng.integers(0, 3))

            plus = mu.copy(); plus[seg, dim] += step
            minus = mu.copy(); minus[seg, dim] -= step
            numeric = (cost(plus) - cost(minus)) / (2 * step)

            assert numeric == pytest.approx(analytic[seg, dim], rel=1e-4, abs=1e-6), (
                "gradient mismatch at segment {} dim {}".format(seg, dim))


class TestWrapper:
    def test_returns_negated_cost_and_flat_gradient(self, problem, reference):
        coords, times, idx_i, idx_j, mu = problem

        loss, grad = cpu_wrapper_chunked(mu.flatten(), coords, times, idx_i, idx_j,
                                         SIGMA, SIGMA_FACTOR)

        assert loss == pytest.approx(-reference[0], rel=1e-10)
        assert grad.shape == (mu.size,)
        np.testing.assert_allclose(grad, -reference[1].flatten(), rtol=1e-9, atol=1e-11)

    def test_both_dispatch_paths_agree(self, problem):
        """The parallel threshold must not change the answer."""
        coords, times, idx_i, idx_j, mu = problem

        serial = cpu_wrapper_chunked(mu.flatten(), coords, times, idx_i, idx_j,
                                     SIGMA, SIGMA_FACTOR, parallel=False)
        par = cpu_wrapper_chunked(mu.flatten(), coords, times, idx_i, idx_j,
                                  SIGMA, SIGMA_FACTOR, parallel=True)

        assert serial[0] == pytest.approx(par[0], rel=1e-10)
        np.testing.assert_allclose(serial[1], par[1], rtol=1e-9, atol=1e-11)

    def test_accepts_the_optimizer_calling_convention(self, problem):
        """The optimizer passes val/deri/chunk_size positionally; None must be fine."""
        coords, times, idx_i, idx_j, mu = problem

        loss, grad = cpu_wrapper_chunked(mu.flatten(), coords, times, idx_i, idx_j,
                                         SIGMA, SIGMA_FACTOR, None, None, int(1e8))

        assert np.isfinite(loss)
        assert np.isfinite(grad).all()

    def test_tolerates_float32_input(self, problem):
        """comet_run_kd hands the CPU backend float32 coordinates."""
        coords, times, idx_i, idx_j, mu = problem

        loss, grad = cpu_wrapper_chunked(mu.flatten(), coords.astype(np.float32),
                                         times.astype(np.int32), idx_i.astype(np.int32),
                                         idx_j.astype(np.int32), SIGMA, SIGMA_FACTOR)

        assert np.isfinite(loss)
        assert np.isfinite(grad).all()

    def test_does_not_mutate_the_drift_estimate(self, problem):
        coords, times, idx_i, idx_j, mu = problem
        flat = mu.flatten()
        original = flat.copy()

        cpu_wrapper_chunked(flat, coords, times, idx_i, idx_j, SIGMA, SIGMA_FACTOR)

        np.testing.assert_array_equal(flat, original)


def test_parallel_threshold_is_a_sane_size():
    # crossover measured at ~200k pairs on an 8-core machine
    assert 50_000 <= PARALLEL_PAIR_THRESHOLD <= 1_000_000
