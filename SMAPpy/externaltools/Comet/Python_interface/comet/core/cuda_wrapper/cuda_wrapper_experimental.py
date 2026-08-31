import math
import numpy as np
from numba import cuda, float64, int32

# ─────────────────────────────────────────────────────────────────────────────
#   TUNE THIS CONSTANTS TO MATCH YOUR ORIGINAL SETTINGS
# ─────────────────────────────────────────────────────────────────────────────

# THREADS_PER_BLOCK should match what you used before (e.g., 128)
THREADS_PER_BLOCK = 128


# ─────────────────────────────────────────────────────────────────────────────
#   REPLACEMENT KERNEL (1D indexing for coord, mu, and deri arrays)
#
#   • “val” accumulation is collapsed into a block‐wide shared‐memory reduction:
#     each block issues exactly one global atomic for val_sum.
#   • Derivatives (d_deri) remain six global atomics per pair (since T is large).
# ─────────────────────────────────────────────────────────────────────────────

@cuda.jit
def cost_function_full_3d_chunked(
    d_locs_time,        # int32[N]
    start_idx,          # int32
    chunk_size,         # int32
    d_idx_i,            # int32[P]
    d_idx_j,            # int32[P]
    d_sigma,            # float64
    d_sigma_factor,     # float64
    d_val_unused,       # float64[P]    (still passed but unused inside)
    d_val_sum,          # float64[1]    (single‐element accumulator for “val”)
    d_deri_flat,        # float64[T*3]  (flattened device array, length = T*3)
    d_locs_coords_flat, # float64[N*3]  (flattened device array, length = N*3)
    mu_flat,            # float64[T*3]  (flattened device array, length = T*3)
    T                   # int32         (number of time‐bins)
):
    """
    Each CUDA block will:
      1) Allocate s_val[THREADS_PER_BLOCK] in shared memory.
      2) Each thread does a grid‐stride loop over (i,j) pairs in [start_idx .. start_idx+chunk_size):
         – compute Gaussian “val” and local derivatives for ti and tj,
         – issue SIX global atomics for derivatives into d_deri_flat[(ti*3)+c] and d_deri_flat[(tj*3)+c],
         – accumulate “val” into a register my_val.
      3) Each thread writes my_val → s_val[threadIdx.x], then we reduce s_val[] → s_val[0].
      4) Thread 0 does ONE global atomic: d_val_sum[0] += s_val[0].
    """

    # Shared‐memory buffer for accumulating “val” per thread
    s_val = cuda.shared.array(shape=THREADS_PER_BLOCK, dtype=float64)

    tx = cuda.threadIdx.x
    bx = cuda.blockIdx.x
    bdim = cuda.blockDim.x
    gdim = cuda.gridDim.x

    # Initialize this thread’s slot in shared “val” array
    s_val[tx] = 0.0
    cuda.syncthreads()

    # Precompute sigma terms once per block (all threads see same d_sigma, d_sigma_factor):
    inv_sigma_scaled = 1.0 / (d_sigma * d_sigma_factor)  # double
    tmp = 2.0 * d_sigma * d_sigma_factor                 # double
    sigma_sq = tmp * tmp                                 # double

    # ─────────────────────────────────────────────────────────────────────────
    # 2) Grid‐stride loop: each thread processes many (i,j) pairs
    # ─────────────────────────────────────────────────────────────────────────
    my_val = float64(0.0)

    pos = bx * bdim + tx
    idx = start_idx + pos
    stride = bdim * gdim

    while pos < chunk_size:
        i = d_idx_i[idx]
        j = d_idx_j[idx]

        ti = d_locs_time[i]
        tj = d_locs_time[j]

        # load coords (flattened): index = 3*i + comp
        xi = d_locs_coords_flat[3 * i + 0]
        yi = d_locs_coords_flat[3 * i + 1]
        zi = d_locs_coords_flat[3 * i + 2]

        xj = d_locs_coords_flat[3 * j + 0]
        yj = d_locs_coords_flat[3 * j + 1]
        zj = d_locs_coords_flat[3 * j + 2]

        # load mu offsets (flattened): index = 3*ti + comp
        mux_ti = mu_flat[3 * ti + 0]
        muy_ti = mu_flat[3 * ti + 1]
        muz_ti = mu_flat[3 * ti + 2]

        mux_tj = mu_flat[3 * tj + 0]
        muy_tj = mu_flat[3 * tj + 1]
        muz_tj = mu_flat[3 * tj + 2]

        # compute displaced differences
        dx = (xi - mux_ti) - (xj - mux_tj)
        dy = (yi - muy_ti) - (yj - muy_tj)
        dz = (zi - muz_ti) - (zj - muz_tj)

        diff_sq = dx * dx + dy * dy + dz * dz
        local_val = inv_sigma_scaled * math.exp(-diff_sq / sigma_sq)
        my_val += local_val

        # derivative factor
        factor = 2.0 * local_val / sigma_sq

        # derivatives for ti
        dx_ti = factor * (  xi - xj + mux_tj - mux_ti )
        dy_ti = factor * (  yi - yj + muy_tj - muy_ti )
        dz_ti = factor * (  zi - zj + muz_tj - muz_ti )

        # derivatives for tj
        dx_tj = factor * (  xj - xi + mux_ti - mux_tj )
        dy_tj = factor * (  yj - yi + muy_ti - muy_tj )
        dz_tj = factor * (  zj - zi + muz_ti - muz_tj )

        # SIX global atomics for derivatives into d_deri_flat[(t*3)+comp]
        base_ti = ti * 3
        cuda.atomic.add(d_deri_flat, base_ti + 0, dx_ti)
        cuda.atomic.add(d_deri_flat, base_ti + 1, dy_ti)
        cuda.atomic.add(d_deri_flat, base_ti + 2, dz_ti)

        base_tj = tj * 3
        cuda.atomic.add(d_deri_flat, base_tj + 0, dx_tj)
        cuda.atomic.add(d_deri_flat, base_tj + 1, dy_tj)
        cuda.atomic.add(d_deri_flat, base_tj + 2, dz_tj)

        # move to next pair
        pos += stride
        idx = start_idx + pos

    cuda.syncthreads()

    # ─────────────────────────────────────────────────────────────────────────
    # 3) Each thread writes my_val → s_val[tx]; then reduce s_val[] → s_val[0]
    # ─────────────────────────────────────────────────────────────────────────
    s_val[tx] = my_val
    cuda.syncthreads()

    offset = bdim >> 1
    while offset > 0:
        if tx < offset:
            s_val[tx] += s_val[tx + offset]
        cuda.syncthreads()
        offset >>= 1

    # ─────────────────────────────────────────────────────────────────────────
    # 4) Thread 0 of each block flushes s_val[0] to d_val_sum[0]
    # ─────────────────────────────────────────────────────────────────────────
    if tx == 0:
        cuda.atomic.add(d_val_sum, 0, s_val[0])








# ─────────────────────────────────────────────────────────────────────────────
#   REPLACEMENT WRAPPER (flattening arrays before passing to kernel)
# ─────────────────────────────────────────────────────────────────────────────

def cuda_wrapper_chunked(
    mu_flat_in,         # 1D np.ndarray, length = T*3, dtype=float64
    d_locs_coords,      # device array float32[N,3]
    d_locs_time,        # device array int32[N]
    d_idx_i,            # device array int32[P]
    d_idx_j,            # device array int32[P]
    d_sigma,            # float64
    d_sigma_factor,     # float64
    d_val,              # device array float64[P]  (still passed but unused by the new kernel)
    d_deri,             # device array float64[T,3] (must be zeroed before call)
    chunk_size          # int
):
    """
    • mu_flat_in: 1D np.ndarray length = T*3 (flattened drift), dtype=float64
    • d_locs_coords: device float32[N,3]
    • d_locs_time:   device int32[N]
    • d_idx_i, d_idx_j: device int32[P]
    • d_sigma, d_sigma_factor: scalars (float64)
    • d_val: device float64[P] (unused inside kernel)
    • d_deri: device float64[T,3] (zeroed before calling this)
    • chunk_size: how many pairs per kernel launch
    """

    P = d_idx_i.size
    # Compute T from length of mu_flat_in
    total_len = mu_flat_in.size
    if total_len % 3 != 0:
        raise ValueError("mu_flat_in length must be a multiple of 3.")
    T = np.int32(total_len // 3)

    # Copy mu_flat_in (already 1D length T*3) to device
    mu_dev = cuda.to_device(mu_flat_in.astype(np.float64))

    # Convert d_locs_coords (float32 N×3 on device) → host → float64 → device (flattened)
    coords_host = d_locs_coords.copy_to_host().astype(np.float64).reshape(-1)
    d_locs_coords_flat = cuda.to_device(coords_host)

    # Flatten d_deri on device: from (T,3) → (T*3,)
    d_deri_flat = d_deri.reshape(d_deri.size)

    # Prepare single‐element accumulator for “val”
    h_val_sum = np.zeros(1, dtype=np.float64)
    d_val_sum = cuda.to_device(h_val_sum)

    threadsperblock = THREADS_PER_BLOCK
    n_chunks = int(math.ceil(P / float(chunk_size)))

    for c in range(n_chunks):
        start_idx = np.int32(c * chunk_size)

        if c == n_chunks - 1:
            this_chunk = np.int32(P - start_idx)
        else:
            this_chunk = np.int32(chunk_size)

        blockspergrid = int(math.ceil(this_chunk / float(threadsperblock)))

        cost_function_full_3d_chunked[blockspergrid, threadsperblock](
            d_locs_time,
            start_idx,
            this_chunk,
            d_idx_i,
            d_idx_j,
            d_sigma,
            d_sigma_factor,
            d_val,              # still passed but unused
            d_val_sum,          # single‐element “val” accumulator
            d_deri_flat,        # flattened device array length = T*3
            d_locs_coords_flat, # flattened device array length = N*3
            mu_dev,             # flattened mu length = T*3
            T
        )

    # Copy back final results
    val_total = d_val_sum.copy_to_host()[0]
    deri_flat_host = d_deri_flat.copy_to_host()   # length = T*3
    deri_host = deri_flat_host.reshape(T, 3)

    # Zero d_deri_flat on device for next iteration
    d_deri_flat[:] = 0.0

    return -val_total, -deri_host.flatten()






