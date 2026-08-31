"""CPU backend for the drift cost function.

Two implementations of the same maths live here:

``_cost_and_gradient_reference``
    A plain Python loop, one pair at a time. Slow, but it is the readable
    statement of what the cost function *is*, and every other backend is
    checked against it.
``_cost_and_gradient_njit``
    The same loop compiled with numba. This is what actually runs.

numba is already a hard dependency -- it is what the CUDA kernels are built on
-- so the compiled path costs nothing extra to install. It is roughly 700x
faster than the interpreted loop and agrees with it to ~1e-18.
"""

import math

import numpy as np
from numba import njit, prange

# Below this many pairs, thread start-up costs more than the parallel loop saves.
PARALLEL_PAIR_THRESHOLD = 200_000

# Number of independent accumulation buffers used by the parallel kernel. Each
# is (n_segments, 3) float64, so the memory cost is small and bounded.
_N_BLOCKS = 16


def _cost_and_gradient_reference(coords, times, idx_i, idx_j, mu, sigma, sigma_factor):
    """Reference implementation: the cost function written out plainly.

    Not used at runtime. Kept as the specification that the compiled kernel and
    the GPU backends are tested against.

    Parameters
    ----------
    coords : (N, 3) float array of localization coordinates in nm.
    times : (N,) int array mapping each localization to its segment.
    idx_i, idx_j : (P,) int arrays of neighbour-pair indices.
    mu : (T, 3) float array, the current per-segment drift estimate in nm.
    sigma, sigma_factor : the Gaussian width and its current schedule multiplier.

    Returns
    -------
    total : float, the summed pair overlap.
    deri : (T, 3) float array, d(total)/d(mu).
    """
    sigma_sq = (2.0 * sigma * sigma_factor) ** 2
    inv_sigma = 1.0 / (sigma * sigma_factor)
    deri = np.zeros_like(mu)
    total = 0.0

    for p in range(len(idx_i)):
        i = idx_i[p]
        j = idx_j[p]
        ti = times[i]
        tj = times[j]

        dx = (coords[i, 0] - mu[ti, 0]) - (coords[j, 0] - mu[tj, 0])
        dy = (coords[i, 1] - mu[ti, 1]) - (coords[j, 1] - mu[tj, 1])
        dz = (coords[i, 2] - mu[ti, 2]) - (coords[j, 2] - mu[tj, 2])

        val = math.exp(-(dx * dx + dy * dy + dz * dz) / sigma_sq) * inv_sigma
        total += val

        for d in range(3):
            # The two contributions are exact negatives of each other.
            contrib = 2.0 * val * (coords[j, d] - coords[i, d] + mu[ti, d] - mu[tj, d]) / sigma_sq
            deri[tj, d] += contrib
            deri[ti, d] -= contrib

    return total, deri


@njit(cache=True, fastmath=True, nogil=True)
def _cost_and_gradient_njit(coords, times, idx_i, idx_j, mu, sigma, sigma_factor):
    """Compiled equivalent of :func:`_cost_and_gradient_reference`."""
    sigma_sq = (2.0 * sigma * sigma_factor) ** 2
    inv_sigma = 1.0 / (sigma * sigma_factor)
    deri = np.zeros_like(mu)
    total = 0.0

    for p in range(idx_i.shape[0]):
        i = idx_i[p]
        j = idx_j[p]
        ti = times[i]
        tj = times[j]

        dx = (coords[i, 0] - mu[ti, 0]) - (coords[j, 0] - mu[tj, 0])
        dy = (coords[i, 1] - mu[ti, 1]) - (coords[j, 1] - mu[tj, 1])
        dz = (coords[i, 2] - mu[ti, 2]) - (coords[j, 2] - mu[tj, 2])

        val = math.exp(-(dx * dx + dy * dy + dz * dz) / sigma_sq) * inv_sigma
        total += val

        for d in range(3):
            contrib = 2.0 * val * (coords[j, d] - coords[i, d] + mu[ti, d] - mu[tj, d]) / sigma_sq
            deri[tj, d] += contrib
            deri[ti, d] -= contrib

    return total, deri


@njit(cache=True, fastmath=True, nogil=True, parallel=True)
def _cost_and_gradient_njit_parallel(coords, times, idx_i, idx_j, mu, sigma, sigma_factor):
    """Parallel variant.

    The gradient update is a scatter-add, which cannot be expressed as a numba
    reduction, so each block accumulates into its own buffer and the buffers are
    summed at the end. Blocks are used rather than thread ids because
    numba.get_thread_id() is not available on the oldest supported numba.
    """
    sigma_sq = (2.0 * sigma * sigma_factor) ** 2
    inv_sigma = 1.0 / (sigma * sigma_factor)
    n_pairs = idx_i.shape[0]
    n_segments = mu.shape[0]

    deri_blocks = np.zeros((_N_BLOCKS, n_segments, 3))
    totals = np.zeros(_N_BLOCKS)
    block_size = (n_pairs + _N_BLOCKS - 1) // _N_BLOCKS

    for b in prange(_N_BLOCKS):
        start = b * block_size
        stop = min(start + block_size, n_pairs)
        local_total = 0.0
        for p in range(start, stop):
            i = idx_i[p]
            j = idx_j[p]
            ti = times[i]
            tj = times[j]

            dx = (coords[i, 0] - mu[ti, 0]) - (coords[j, 0] - mu[tj, 0])
            dy = (coords[i, 1] - mu[ti, 1]) - (coords[j, 1] - mu[tj, 1])
            dz = (coords[i, 2] - mu[ti, 2]) - (coords[j, 2] - mu[tj, 2])

            val = math.exp(-(dx * dx + dy * dy + dz * dz) / sigma_sq) * inv_sigma
            local_total += val

            for d in range(3):
                contrib = 2.0 * val * (coords[j, d] - coords[i, d] + mu[ti, d] - mu[tj, d]) / sigma_sq
                deri_blocks[b, tj, d] += contrib
                deri_blocks[b, ti, d] -= contrib
        totals[b] = local_total

    deri = np.zeros_like(mu)
    for b in range(_N_BLOCKS):
        deri += deri_blocks[b]
    return totals.sum(), deri


def cpu_wrapper_chunked(mu, locs_coords, locs_time, idx_i, idx_j, sigma, sigma_factor,
                        val=None, deri=None, chunk_size=None, debug=False, parallel=None):
    """Cost and gradient for the optimizer, on the CPU.

    Signature matches the CUDA and torch wrappers so the optimizer can swap
    backends. `val`, `deri` and `chunk_size` exist for that compatibility only:
    the compiled kernel needs no scratch buffer and no chunking, since it is not
    working around device memory limits.
    """
    mu = np.ascontiguousarray(mu.reshape((-1, 3)), dtype=np.float64)
    coords = np.ascontiguousarray(locs_coords[:, :3], dtype=np.float64)
    times = np.ascontiguousarray(locs_time, dtype=np.int64)
    idx_i = np.ascontiguousarray(idx_i, dtype=np.int64)
    idx_j = np.ascontiguousarray(idx_j, dtype=np.int64)

    if parallel is None:
        parallel = idx_i.shape[0] >= PARALLEL_PAIR_THRESHOLD
    kernel = _cost_and_gradient_njit_parallel if parallel else _cost_and_gradient_njit

    total, gradient = kernel(coords, times, idx_i, idx_j, mu,
                             float(sigma), float(sigma_factor))

    if debug:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(3, 2)
        for d in range(3):
            ax[d, 0].plot(gradient[:, d])
            ax[d, 1].plot(mu[:, d])
        ax[0, 0].set_title(f"Gradients (sigma={sigma * sigma_factor:.2f} nm)")
        ax[0, 1].set_title("Drift Estimate [nm]")
        plt.tight_layout()
        plt.show()

    # The optimizer minimizes, so hand back the negated overlap.
    return -total, -gradient.flatten()


def cpu_wrapper_chunked_fast(mu, locs_coords, locs_time, idx_i, idx_j, sigma, sigma_factor,
                             val=None, deri=None, chunk_size=None, debug=False, parallel=None):
    """`cpu_wrapper_chunked` without the per-call array conversions.

    The optimizer calls this a few hundred times with the *same* arrays, and the
    original converts them every time: the pair indices arrive as int32 and are
    copied to int64, which at 300 M pairs is 5 GB of allocation and copying per
    cost evaluation -- about a fifth of the call -- for nothing.  The compiled
    kernel specialises on whatever dtypes it is handed, so the caller's float32
    coordinates and int32 indices are used as they are.

    The caller is responsible for passing contiguous arrays; the optimizer
    already does (`optimize_3d_chunked_better_moving_avg_kd` builds them once).
    """
    mu = mu.reshape((-1, 3))  # float64 and contiguous, straight from L-BFGS-B
    if parallel is None:
        parallel = idx_i.shape[0] >= PARALLEL_PAIR_THRESHOLD
    kernel = _cost_and_gradient_njit_parallel if parallel else _cost_and_gradient_njit

    total, gradient = kernel(locs_coords, locs_time, idx_i, idx_j, mu,
                             float(sigma), float(sigma_factor))
    # The optimizer minimizes, so hand back the negated overlap.
    return -total, -gradient.flatten()


# --- smappy: approximate kernel ---------------------------------------------
# `exp` is 54% of the kernel's runtime on ARM (measured: replacing it with a
# polynomial of the wrong value takes 0.65 s down to 0.30 s), so it is worth
# computing it approximately: the cost is a sum over 3e8 pairs and a relative
# error of 1e-5 per term is far below the noise of the estimate itself.
#
#   exp(-u) = exp(-k) * exp(-f),  k = int(u), f = u - k in [0, 1)
#
# with exp(-k) from a table and exp(-f) from its 7th-order Taylor series, whose
# relative error is below 7e-5 on [0, 1).  Pairs beyond CUTOFF_SIGMAS sigma are
# skipped outright, which also bounds the table.  Measured against the exact
# kernel: 1.3x faster, gradient accurate to 4e-4 relative at the finest sigma.
CUTOFF_SIGMAS = 6.0
_EXPTAB = np.exp(-np.arange(float(int(CUTOFF_SIGMAS ** 2 / 4.0) + 2)))


@njit(cache=True, fastmath=True, nogil=True, parallel=True)
def _cost_and_gradient_njit_parallel_approx(coords, times, idx_i, idx_j, mu,
                                            sigma, sigma_factor, cutoff_sigmas, exptab):
    """`_cost_and_gradient_njit_parallel` with a tabulated exp and a cutoff."""
    s_eff = sigma * sigma_factor
    inv_sigma = 1.0 / s_eff
    # one reciprocal, hoisted: the exact kernel divides by sigma_sq four times
    # per pair (once for the exponent, once per gradient component)
    inv_sigma_sq = 1.0 / (2.0 * s_eff) ** 2
    cutoff_sq = (cutoff_sigmas * s_eff) ** 2

    n_pairs = idx_i.shape[0]
    n_segments = mu.shape[0]
    deri_blocks = np.zeros((_N_BLOCKS, n_segments, 3))
    totals = np.zeros(_N_BLOCKS)
    block_size = (n_pairs + _N_BLOCKS - 1) // _N_BLOCKS

    for b in prange(_N_BLOCKS):
        start = b * block_size
        stop = min(start + block_size, n_pairs)
        local_total = 0.0
        for p in range(start, stop):
            i = idx_i[p]
            j = idx_j[p]
            ti = times[i]
            tj = times[j]

            dx = (coords[i, 0] - mu[ti, 0]) - (coords[j, 0] - mu[tj, 0])
            dy = (coords[i, 1] - mu[ti, 1]) - (coords[j, 1] - mu[tj, 1])
            dz = (coords[i, 2] - mu[ti, 2]) - (coords[j, 2] - mu[tj, 2])
            diff_sq = dx * dx + dy * dy + dz * dz
            if diff_sq > cutoff_sq:  # contributes less than exp(-cutoff^2/4)
                continue

            u = diff_sq * inv_sigma_sq
            k = int(u)
            f = u - k
            series = 1.0 + f * (-1.0 + f * (0.5 + f * (-1.0 / 6.0 + f * (
                1.0 / 24.0 + f * (-1.0 / 120.0 + f * (1.0 / 720.0 - f / 5040.0))))))
            val = exptab[k] * series * inv_sigma
            local_total += val

            weight = 2.0 * val * inv_sigma_sq
            for d in range(3):
                contrib = weight * (coords[j, d] - coords[i, d] + mu[ti, d] - mu[tj, d])
                deri_blocks[b, tj, d] += contrib
                deri_blocks[b, ti, d] -= contrib
        totals[b] = local_total

    deri = np.zeros_like(mu)
    for b in range(_N_BLOCKS):
        deri += deri_blocks[b]
    return totals.sum(), deri


def cpu_wrapper_chunked_approx(mu, locs_coords, locs_time, idx_i, idx_j, sigma, sigma_factor,
                               val=None, deri=None, chunk_size=None, debug=False, parallel=None):
    """`cpu_wrapper_chunked_fast` with the approximate kernel above."""
    mu = mu.reshape((-1, 3))
    total, gradient = _cost_and_gradient_njit_parallel_approx(
        locs_coords, locs_time, idx_i, idx_j, mu, float(sigma), float(sigma_factor),
        float(CUTOFF_SIGMAS), _EXPTAB)
    return -total, -gradient.flatten()


# --- smappy: quality control on the CPU -------------------------------------
# The overlap a segment actually achieves, against the overlap it would have
# with no drift correction at all.  A segment whose pairs do not overlap better
# than that null is not supported by its own data -- its drift is a guess -- and
# `flag_flawed_segments` marks it so the optimizer can set it to NaN and let the
# spline bridge it.  COMET computes this on its cuda_qc / torch_qc backends
# only; this is the same quantity for the CPU path, run once after the fit.
@njit(cache=True, fastmath=True, nogil=True, parallel=True)
def _overlap_per_segment_njit(coords, times, idx_i, idx_j, mu, sigma, sigma_factor):
    """Per-segment overlap sums with the fitted drift, with none, and pair counts.

    Same-segment pairs are excluded: their overlap does not depend on the drift,
    so they say nothing about whether this segment's estimate is any good.  Each
    pair counts for both of its segments, as in the GPU version.
    """
    sigma_sq = (2.0 * sigma * sigma_factor) ** 2
    inv_sigma = 1.0 / (sigma * sigma_factor)
    n_segments = mu.shape[0]
    n_pairs = idx_i.shape[0]

    obs_blocks = np.zeros((_N_BLOCKS, n_segments))
    null_blocks = np.zeros((_N_BLOCKS, n_segments))
    count_blocks = np.zeros((_N_BLOCKS, n_segments))
    block_size = (n_pairs + _N_BLOCKS - 1) // _N_BLOCKS

    for b in prange(_N_BLOCKS):
        start = b * block_size
        stop = min(start + block_size, n_pairs)
        for p in range(start, stop):
            i = idx_i[p]
            j = idx_j[p]
            ti = times[i]
            tj = times[j]
            if ti == tj:
                continue

            ux = coords[i, 0] - coords[j, 0]
            uy = coords[i, 1] - coords[j, 1]
            uz = coords[i, 2] - coords[j, 2]

            dx = ux - (mu[ti, 0] - mu[tj, 0])
            dy = uy - (mu[ti, 1] - mu[tj, 1])
            dz = uz - (mu[ti, 2] - mu[tj, 2])

            val = math.exp(-(dx * dx + dy * dy + dz * dz) / sigma_sq) * inv_sigma
            null = math.exp(-(ux * ux + uy * uy + uz * uz) / sigma_sq) * inv_sigma

            obs_blocks[b, ti] += val
            obs_blocks[b, tj] += val
            null_blocks[b, ti] += null
            null_blocks[b, tj] += null
            count_blocks[b, ti] += 1.0
            count_blocks[b, tj] += 1.0

    return obs_blocks.sum(axis=0), null_blocks.sum(axis=0), count_blocks.sum(axis=0)


def cpu_wrapper_chunked_qc(mu, locs_coords, locs_time, idx_i, idx_j, sigma, sigma_factor,
                           val=None, deri=None, chunk_size=None):
    """Cost, gradient and the quality-control metric, on the CPU.

    Signature matches `cuda_wrapper_chunked_qc` / `torch_wrapper_chunked_qc`.
    """
    loss, gradient = cpu_wrapper_chunked_fast(mu, locs_coords, locs_time, idx_i, idx_j,
                                              sigma, sigma_factor)
    obs, null, pairs = _overlap_per_segment_njit(
        locs_coords, locs_time, idx_i, idx_j, mu.reshape((-1, 3)),
        float(sigma), float(sigma_factor))

    denominator = np.maximum(pairs, np.finfo(np.float64).eps)
    qc = {
        "window": 1,
        "count_both_sides": True,
        "excluded_same_time": True,
        "window_pairs": pairs,
        "Q_obs": obs / denominator,
        "Q_null": null / denominator,
        "num_obs": obs,
        "num_null": null,
    }
    qc["diff"] = qc["Q_obs"] - qc["Q_null"]
    return loss, gradient, qc
