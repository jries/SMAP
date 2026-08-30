import math
import numpy as np
from numba import cuda


@cuda.jit
def cost_function_full_3d_chunked_qc(
        d_locs_time,  # (N,) int32
        start_idx,  # int
        chunk_size,  # int
        d_idx_i, d_idx_j,  # (P,) int32
        inv_sigma,  # float64
        sigma_sq,  # float64
        d_locs_coords,  # (N,3) float64
        mu,  # (T,3) float64
        d_deri,  # (T,3) float64 (atomic add)
        d_per_overlap_obs,  # (T,)  float64 (atomic add)
        d_per_overlap_null,  # (T,)  float64 (atomic add)
        d_per_pairs,  # (T,)  float64 (atomic add)
        exclude_same_time,  # int32 (0/1)
        count_both_sides,  # int32 (0/1)
        d_val_sum  # (1,)  float64 (atomic add) for loss
):
    """
    Fused cost + QC per-frame accumulators.
    - observed overlap uses current mu
    - null   overlap uses mu == 0
    - pairs  counts how many pairs "touch" each frame (ti and optionally tj)
    """
    tx = cuda.threadIdx.x
    ty = cuda.blockIdx.x
    bw = cuda.blockDim.x
    pos = tx + ty * bw

    if pos >= chunk_size:
        return

    pair_idx = pos + start_idx
    i = d_idx_i[pair_idx]
    j = d_idx_j[pair_idx]

    ti = d_locs_time[i]
    tj = d_locs_time[j]

    # optionally exclude same-time pairs
    if exclude_same_time and (ti == tj):
        return

    # prefetch coords
    xi0 = d_locs_coords[i, 0];
    xi1 = d_locs_coords[i, 1];
    xi2 = d_locs_coords[i, 2]
    xj0 = d_locs_coords[j, 0];
    xj1 = d_locs_coords[j, 1];
    xj2 = d_locs_coords[j, 2]

    # Δx_ij
    dd0 = xi0 - xj0
    dd1 = xi1 - xj1
    dd2 = xi2 - xj2

    # === observed (with mu)
    mu_ti0 = mu[ti, 0];
    mu_ti1 = mu[ti, 1];
    mu_ti2 = mu[ti, 2]
    mu_tj0 = mu[tj, 0];
    mu_tj1 = mu[tj, 1];
    mu_tj2 = mu[tj, 2]
    # d_obs = (Δx_ij) - (mu_ti - mu_tj)
    do0 = dd0 - (mu_ti0 - mu_tj0)
    do1 = dd1 - (mu_ti1 - mu_tj1)
    do2 = dd2 - (mu_ti2 - mu_tj2)
    diff_sq_obs = do0 * do0 + do1 * do1 + do2 * do2
    val_obs = inv_sigma * math.exp(-diff_sq_obs / sigma_sq)

    # === null (mu == 0): d0 = Δx_ij
    diff_sq_null = dd0 * dd0 + dd1 * dd1 + dd2 * dd2
    val_null = inv_sigma * math.exp(-diff_sq_null / sigma_sq)

    # ---- accumulate loss (observed only) and gradient (same as your code)
    cuda.atomic.add(d_val_sum, 0, val_obs)

    # grad wrt mu[tj] and mu[ti]
    # contrib_tj = 2*val_obs*(xj - xi + mu_ti - mu_tj)/sigma_sq
    c0 = 2.0 * val_obs * (xj0 - xi0 + mu_ti0 - mu_tj0) / sigma_sq
    c1 = 2.0 * val_obs * (xj1 - xi1 + mu_ti1 - mu_tj1) / sigma_sq
    c2 = 2.0 * val_obs * (xj2 - xi2 + mu_ti2 - mu_tj2) / sigma_sq
    cuda.atomic.add(d_deri, (tj, 0), c0)
    cuda.atomic.add(d_deri, (tj, 1), c1)
    cuda.atomic.add(d_deri, (tj, 2), c2)

    # contrib_ti = 2*val_obs*(xi - xj + mu_tj - mu_ti)/sigma_sq
    c0 = 2.0 * val_obs * (xi0 - xj0 + mu_tj0 - mu_ti0) / sigma_sq
    c1 = 2.0 * val_obs * (xi1 - xj1 + mu_tj1 - mu_ti1) / sigma_sq
    c2 = 2.0 * val_obs * (xi2 - xj2 + mu_tj2 - mu_ti2) / sigma_sq
    cuda.atomic.add(d_deri, (ti, 0), c0)
    cuda.atomic.add(d_deri, (ti, 1), c1)
    cuda.atomic.add(d_deri, (ti, 2), c2)

    # ---- QC numerators: per-frame overlap (observed + null)
    cuda.atomic.add(d_per_overlap_obs, ti, val_obs)
    cuda.atomic.add(d_per_overlap_null, ti, val_null)
    if count_both_sides:
        cuda.atomic.add(d_per_overlap_obs, tj, val_obs)
        cuda.atomic.add(d_per_overlap_null, tj, val_null)

    # ---- QC denominator: per-frame pair counts
    cuda.atomic.add(d_per_pairs, ti, 1.0)
    if count_both_sides:
        cuda.atomic.add(d_per_pairs, tj, 1.0)


def _moving_window_sum_1d(x: np.ndarray, w: int) -> np.ndarray:
    """Centered box filter with zero padding (like PyTorch conv1d padding='same')."""
    if w is None or w <= 1:
        return x.astype(np.float64, copy=True)
    ker = np.ones(int(w), dtype=np.float64)
    return np.convolve(x.astype(np.float64), ker, mode="same")


def cuda_wrapper_chunked_qc(
        mu, d_locs_coords, d_locs_time, d_idx_i, d_idx_j, d_sigma, d_sigma_factor, d_val, d_deri, chunk_size):
    """
    Launches the fused QC kernel over chunks and computes windowed, pairs-normalized metrics
    against a zero-drift baseline.

    Returns:
      loss: float
      grad_flat: np.ndarray (T*3,)
      qc: dict with Q_obs, Q_null, diff, lift, window_pairs etc.
    """
    # ---- prepare scalars (float64 for stability)
    sigma = float(d_sigma)
    sigma_f = float(d_sigma_factor)
    inv_sigma = 1.0 / (sigma * sigma_f)
    sigma_sq = (2.0 * sigma * sigma_f) ** 2

    window = 1
    exclude_same_time = True
    count_both_sides = True

    # ---- device buffers
    T = int(mu.size // 3)
    mu_dev = cuda.to_device(np.asarray(mu, dtype=np.float64).reshape(T, 3))

    d_deri = cuda.to_device(np.zeros((T, 3), dtype=np.float64))
    d_val_sum = cuda.to_device(np.zeros(1, dtype=np.float64))

    per_overlap_obs_dev = cuda.to_device(np.zeros(T, dtype=np.float64))
    per_overlap_null_dev = cuda.to_device(np.zeros(T, dtype=np.float64))
    per_pairs_dev = cuda.to_device(np.zeros(T, dtype=np.float64))

    # ---- grid
    threadsperblock = 128
    n_pairs = int(d_idx_i.size)
    n_chunks = int(np.ceil(n_pairs / chunk_size))

    # ---- loop
    for c in range(n_chunks):
        n_rem = min(chunk_size, n_pairs - c * chunk_size)
        if n_rem <= 0:
            break
        idc_start = c * chunk_size
        blockspergrid = (n_rem + (threadsperblock - 1)) // threadsperblock

        cost_function_full_3d_chunked_qc[blockspergrid, threadsperblock](
            d_locs_time,
            idc_start, n_rem,
            d_idx_i, d_idx_j,
            inv_sigma, sigma_sq,
            d_locs_coords, mu_dev,
            d_deri,
            per_overlap_obs_dev, per_overlap_null_dev,
            per_pairs_dev,
            np.int32(1 if exclude_same_time else 0),
            np.int32(1 if count_both_sides else 0),
            d_val_sum
        )

    # ---- fetch results
    loss = -float(d_val_sum.copy_to_host()[0])
    deri = d_deri.copy_to_host();
    d_deri[:] = 0
    grad_flat = (-deri).reshape(-1)

    per_overlap_obs = per_overlap_obs_dev.copy_to_host()
    per_overlap_null = per_overlap_null_dev.copy_to_host()
    per_pairs = per_pairs_dev.copy_to_host()

    # ---- windowed sums (CPU)
    w = 1 if (window is None or window <= 1) else int(window)
    w_obs = _moving_window_sum_1d(per_overlap_obs, w)
    w_null = _moving_window_sum_1d(per_overlap_null, w)
    w_pairs = _moving_window_sum_1d(per_pairs, w)

    # ---- normalized curves
    eps = np.finfo(np.float64).eps
    denom = np.maximum(w_pairs, eps)
    Q_obs = w_obs / denom
    Q_null = w_null / denom
    diff = Q_obs - Q_null
    lift = (Q_obs / np.maximum(Q_null, eps)) - 1.0

    qc = {
        "window": w,
        "excluded_same_time": bool(exclude_same_time),
        "count_both_sides": bool(count_both_sides),
        "window_pairs": w_pairs,
        "Q_obs": Q_obs,
        "Q_null": Q_null,
        "diff": diff,
        "lift": lift,
        "num_obs": w_obs,
        "num_null": w_null,
    }
    return loss, grad_flat, qc
