import math
import numpy as np
from numba import cuda
import matplotlib.pyplot as plt

debug_mode = False  # Set True for visual debugging


# Numba cuda code of cost-function --> quantized overlap of pairs of localizations, gaussian mixture model
@cuda.jit
def cost_function_full_3d_chunked(d_locs_time, start_idx, chunk_size, d_idx_i, d_idx_j, d_sigma, d_sigma_factor, d_val,
                                  d_val_sum, d_deri, d_locs_coords, mu):
    """Compute negative-overlap cost and gradient for 3D localizations (GPU/CPU backend; internal use)."""
    tx = cuda.threadIdx.x
    ty = cuda.blockIdx.x
    bw = cuda.blockDim.x
    pos = tx + ty * bw

    if pos < chunk_size:
        i = d_idx_i[pos + start_idx]
        j = d_idx_j[pos + start_idx]

        ti = d_locs_time[i]
        tj = d_locs_time[j]

        dx = (d_locs_coords[i, 0] - mu[ti, 0]) - (d_locs_coords[j, 0] - mu[tj, 0])
        dy = (d_locs_coords[i, 1] - mu[ti, 1]) - (d_locs_coords[j, 1] - mu[tj, 1])
        dz = (d_locs_coords[i, 2] - mu[ti, 2]) - (d_locs_coords[j, 2] - mu[tj, 2])
        sigma_sq = (2 * d_sigma * d_sigma_factor) ** 2

        diff_sq = dx * dx + dy * dy + dz * dz
        val = 1 / (d_sigma * d_sigma_factor) * math.exp(-diff_sq / sigma_sq)
        d_val[pos] = val

        # Update derivatives
        cuda.atomic.add(d_deri, (tj, 0), 2 * val * (d_locs_coords[j, 0] - d_locs_coords[i, 0] + mu[ti, 0] - mu[tj, 0]) / sigma_sq)
        cuda.atomic.add(d_deri, (tj, 1), 2 * val * (d_locs_coords[j, 1] - d_locs_coords[i, 1] + mu[ti, 1] - mu[tj, 1]) / sigma_sq)
        cuda.atomic.add(d_deri, (tj, 2), 2 * val * (d_locs_coords[j, 2] - d_locs_coords[i, 2] + mu[ti, 2] - mu[tj, 2]) / sigma_sq)

        cuda.atomic.add(d_deri, (ti, 0), 2 * val * (d_locs_coords[i, 0] - d_locs_coords[j, 0] + mu[tj, 0] - mu[ti, 0]) / sigma_sq)
        cuda.atomic.add(d_deri, (ti, 1), 2 * val * (d_locs_coords[i, 1] - d_locs_coords[j, 1] + mu[tj, 1] - mu[ti, 1]) / sigma_sq)
        cuda.atomic.add(d_deri, (ti, 2), 2 * val * (d_locs_coords[i, 2] - d_locs_coords[j, 2] + mu[tj, 2] - mu[ti, 2]) / sigma_sq)

        cuda.atomic.add(d_val_sum, 0, val)
        d_val[pos] = 0


# Interface between the Python code and the CUDA kernel, mainly for chunking the data to avoid memory issues
def cuda_wrapper_chunked(mu, d_locs_coords, d_locs_time, d_idx_i, d_idx_j, d_sigma, d_sigma_factor, d_val, d_deri,
                         chunk_size):
    d_val_sum = cuda.to_device(np.zeros(1, dtype=np.float64))
    mu_dev = cuda.to_device(np.asarray(mu.reshape(int(mu.size / 3), 3), dtype=np.float64))

    n_chunks = int(np.ceil(d_idx_i.size / chunk_size))
    threadsperblock = 128

    for i in range(n_chunks - 1):
        idc_start = i*chunk_size
        blockspergrid = (chunk_size + (threadsperblock - 1)) // threadsperblock
        cost_function_full_3d_chunked[blockspergrid, threadsperblock](
            d_locs_time, idc_start, chunk_size, d_idx_i, d_idx_j,
            d_sigma, d_sigma_factor, d_val, d_val_sum, d_deri, d_locs_coords, mu_dev
        )

    # Final chunk
    n_remaining = d_idx_i.size - (n_chunks - 1) * chunk_size
    idc_start = (n_chunks - 1) * chunk_size
    blockspergrid = (n_remaining + (threadsperblock - 1)) // threadsperblock
    cost_function_full_3d_chunked[blockspergrid, threadsperblock](
        d_locs_time, idc_start, n_remaining, d_idx_i, d_idx_j,
        d_sigma, d_sigma_factor, d_val, d_val_sum, d_deri, d_locs_coords, mu_dev
    )
    # d_val_sum is a device-side accumulator that is never reset between chunks,
    # so it already holds the total over every chunk. Adding the intermediate
    # copies as well would weight chunk k by (n_chunks - k + 1), handing L-BFGS-B
    # a cost that disagrees with the gradient.
    val_total = d_val_sum.copy_to_host()
    deri = d_deri.copy_to_host()
    d_deri[:] = 0

    if debug_mode:
        mu_host = mu_dev.copy_to_host()
        fig, ax = plt.subplots(3, 2)
        ax[0, 0].plot(deri[:, 0])
        ax[1, 0].plot(deri[:, 1])
        ax[2, 0].plot(deri[:, 2])
        ax[0, 0].set_title(f"Derivatives (σ = {d_sigma * d_sigma_factor} nm)")
        ax[0, 1].plot(mu_host[:, 0])
        ax[1, 1].plot(mu_host[:, 1])
        ax[2, 1].plot(mu_host[:, 2])
        ax[0, 1].set_title("Drift μ")
        plt.tight_layout()
        plt.show()

    return -np.nansum(val_total), -deri.flatten()
