import numpy as np
import torch


def cost_function_full_3d_chunked_torch(
    mu,
    locs_coords,
    locs_time,
    idx_i, idx_j,
    sigma, sigma_f,
    chunk_size=100_000_000,
):
    device, dtype = mu.device, mu.dtype
    inv_sigma = 1.0 / (sigma * sigma_f)
    sigma_sq = (2.0 * sigma * sigma_f) ** 2

    T = mu.shape[0]
    deri = torch.zeros((T, 3), dtype=dtype, device=device)
    total_val = torch.zeros((), dtype=dtype, device=device)

    P = idx_i.numel()
    for s in range(0, P, chunk_size):
        e = min(P, s + chunk_size)
        i = idx_i[s:e]
        j = idx_j[s:e]
        ti = locs_time[i]
        tj = locs_time[j]

        ri = locs_coords[i] - mu[ti]  # (K,3)
        rj = locs_coords[j] - mu[tj]  # (K,3)
        d = ri - rj
        diff_sq = (d * d).sum(dim=1)  # (K,)
        val = inv_sigma * torch.exp(-diff_sq / sigma_sq)

        total_val = total_val + val.sum()

        contrib_tj = 2.0 * val[:, None] * (locs_coords[j] - locs_coords[i] + mu[ti] - mu[tj]) / sigma_sq
        contrib_ti = 2.0 * val[:, None] * (locs_coords[i] - locs_coords[j] + mu[tj] - mu[ti]) / sigma_sq

        deri.index_add_(0, tj, contrib_tj)
        deri.index_add_(0, ti, contrib_ti)

    loss = -total_val
    grad_flat = (-deri).reshape(-1).detach().cpu().numpy()
    return float(loss.item()), grad_flat


def torch_wrapper_chunked(mu, d_locs_coords, d_locs_time,
                          d_idx_i, d_idx_j, d_sigma, sigma_factor, d_val, _,
                          chunk_size):
    device = d_val  # see drift optimizer: d_val is not used in pytorch implementation, but abused as device indicator
    mu_np = np.reshape(mu, (-1, 3))  # (T,3)
    d_mu = torch.as_tensor(mu_np, dtype=torch.float32, device=device)
    d_sigma_factor = torch.as_tensor(sigma_factor, dtype=torch.float32, device=device)
    loss, grad = cost_function_full_3d_chunked_torch(
        d_mu, d_locs_coords, d_locs_time, d_idx_i, d_idx_j, d_sigma, d_sigma_factor, chunk_size
    )
    return loss, grad.astype(np.float64)
