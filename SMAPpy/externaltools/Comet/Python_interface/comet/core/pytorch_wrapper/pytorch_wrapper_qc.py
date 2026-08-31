import torch.nn.functional as F
import torch
import numpy as np
import matplotlib.pyplot as plt

from comet.core.qc_utils import plot_q_with_baseline


@torch.no_grad()
def cost_with_overlap_window_and_zero_null_torch(
    mu, locs_coords, locs_time, idx_i, idx_j, sigma, sigma_factor,
    chunk_size=10_000_000,
    window=1,                   # odd int for centered window; 1/None/1 => no window
    exclude_same_time=True,
    count_both_sides=True,       # numerator & denominator count rule
):
    """
    Loss/grad (as in your Torch cost) + QC metric:
      windowed overlap per time normalized by windowed *pair count*,
    plus a zero-drift (μ=0) null baseline.
    """
    device, dtype = mu.device, mu.dtype
    mu = mu.reshape(-1, 3)
    T = mu.shape[0]

    sigma = torch.as_tensor(sigma, dtype=dtype, device=device)
    sigma_f = torch.as_tensor(sigma_factor, dtype=dtype, device=device)
    inv_sigma = 1.0 / (sigma * sigma_f)
    sigma_sq = (2.0 * sigma * sigma_f) ** 2
    eps = torch.finfo(dtype).eps

    # === loss/grad accumulators (observed)
    deri = torch.zeros((T, 3), dtype=dtype, device=device)
    total_val = torch.zeros((), dtype=dtype, device=device)

    # === per-frame accumulators
    per_frame_overlap_obs  = torch.zeros((T,), dtype=dtype, device=device)
    per_frame_overlap_null = torch.zeros((T,), dtype=dtype, device=device)
    per_frame_pairs        = torch.zeros((T,), dtype=dtype, device=device)  # float for conv1d

    P = idx_i.numel()
    for s in range(0, P, chunk_size):
        e = min(P, s + chunk_size)
        i  = idx_i[s:e]
        j  = idx_j[s:e]
        ti = locs_time[i]
        tj = locs_time[j]

        mask = (ti != tj) if exclude_same_time else torch.ones_like(ti, dtype=torch.bool)
        if not torch.any(mask):
            continue
        i  = i[mask];  j  = j[mask]
        ti = ti[mask]; tj = tj[mask]

        # Precompute pair deltas once: Δx_ij
        delta_ij = locs_coords[i] - locs_coords[j]  # (K,3)

        # === observed μ: d = (Δx_ij) - (μ_ti - μ_tj)
        mu_diff = mu[ti] - mu[tj]
        d_obs   = delta_ij - mu_diff
        diff_sq_obs = (d_obs * d_obs).sum(dim=1)
        val_obs = inv_sigma * torch.exp(-diff_sq_obs / sigma_sq)

        # ---- loss/grad (unchanged)
        total_val += val_obs.sum()
        contrib_tj = 2.0 * val_obs[:, None] * (locs_coords[j] - locs_coords[i] + mu[ti] - mu[tj]) / sigma_sq
        contrib_ti = 2.0 * val_obs[:, None] * (locs_coords[i] - locs_coords[j] + mu[tj] - mu[ti]) / sigma_sq
        deri.index_add_(0, tj, contrib_tj)
        deri.index_add_(0, ti, contrib_ti)

        # accumulate observed per-frame overlap
        per_frame_overlap_obs.index_add_(0, ti, val_obs)
        if count_both_sides:
            per_frame_overlap_obs.index_add_(0, tj, val_obs)

        # === zero-drift null (μ ≡ 0): d0 = Δx_ij
        diff_sq_null = (delta_ij * delta_ij).sum(dim=1)
        val_null = inv_sigma * torch.exp(-diff_sq_null / sigma_sq)

        per_frame_overlap_null.index_add_(0, ti, val_null)
        if count_both_sides:
            per_frame_overlap_null.index_add_(0, tj, val_null)

        # per-frame pair counts (same counting rule)
        ones = torch.ones_like(val_obs, dtype=dtype, device=device)
        per_frame_pairs.index_add_(0, ti, ones)
        if count_both_sides:
            per_frame_pairs.index_add_(0, tj, ones)

    # finalize loss/grad
    loss = -total_val
    grad_flat = (-deri).reshape(-1).detach().cpu().numpy()

    # windowed sums via 1D box filter
    if window is None or window <= 1:
        window = 1
        w_obs   = per_frame_overlap_obs
        w_null  = per_frame_overlap_null
        w_pairs = per_frame_pairs
    else:
        wker = torch.ones((1, 1, int(window)), dtype=dtype, device=device)
        pad  = int(window) // 2
        w_obs   = F.conv1d(per_frame_overlap_obs.view(1,1,-1),  wker, padding=pad).view(-1)
        w_null  = F.conv1d(per_frame_overlap_null.view(1,1,-1), wker, padding=pad).view(-1)
        w_pairs = F.conv1d(per_frame_pairs.view(1,1,-1),        wker, padding=pad).view(-1)

    # denominator (pairs-normalized)
    denom = torch.clamp(w_pairs, min=eps)

    Q_obs  = w_obs  / denom
    Q_null = w_null / denom

    # effect sizes (no null variance here; use lift & absolute diff)
    diff = Q_obs - Q_null
    lift = (Q_obs / torch.clamp(Q_null, min=eps)) - 1.0  # relative improvement

    qc = {
        "window": int(window),
        "count_both_sides": bool(count_both_sides),
        "excluded_same_time": bool(exclude_same_time),
        "window_pairs": w_pairs.detach().cpu().numpy(),
        "Q_obs":  Q_obs.detach().cpu().numpy(),
        "Q_null": Q_null.detach().cpu().numpy(),
        "diff":   diff.detach().cpu().numpy(),
        "lift":   lift.detach().cpu().numpy(),

        # raw numerators (can be useful to inspect)
        "num_obs":  w_obs.detach().cpu().numpy(),
        "num_null": w_null.detach().cpu().numpy(),
    }
    return float(loss.item()), grad_flat, qc


def torch_wrapper_chunked_qc(
    mu, d_locs_coords, d_locs_time, d_idx_i, d_idx_j, d_sigma, d_sigma_factor, d_val, d_deri, chunk_size
):
    device = d_locs_coords.device
    mu_np = np.reshape(mu, (-1, 3))
    d_mu = torch.as_tensor(mu_np, dtype=d_locs_coords.dtype, device=device)

    loss, grad, qc = cost_with_overlap_window_and_zero_null_torch(
        d_mu, d_locs_coords, d_locs_time, d_idx_i, d_idx_j,
        d_sigma, d_sigma_factor,
        chunk_size=chunk_size,
        window=1,
        exclude_same_time=True,
        count_both_sides=True,
    )

    return loss, grad.astype(np.float64), qc