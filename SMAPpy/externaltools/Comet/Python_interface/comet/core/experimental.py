# temporal_refined_drift.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.ndimage import convolve
from numba import cuda

# Import your CUDA and segmentation logic
from comet.core.segmenter import segmentation_wrapper
from comet.core.interpolation import interpolate_drift
from comet.core.pair_indices import pair_indices_kdtree
from comet.core.cuda_wrapper import cuda_wrapper_chunked
from comet.core.io_utils import load_simulation_dataset_and_gt_drift


def compute_gradient_difference(drift_a, drift_b):
    """Mean squared difference between two interpolated drift trajectories."""
    return np.mean(np.linalg.norm(drift_a - drift_b, axis=1))


def optimize_drift(
    locs_nm,
    loc_frames,
    n_segments,
    sigma_nm,
    drift_init=None,
    boxcar_width=3,
    drift_max_nm=300,
    drift_max_bound_factor=2,
    d_coords=None,
    d_times=None,
    d_idx_i=None,
    d_idx_j=None,
    chunk_size=int(1e8),
):
    drift_est = drift_init if drift_init is not None else np.zeros(n_segments * 3)
    bounds = [
        (-drift_max_nm * drift_max_bound_factor, drift_max_nm * drift_max_bound_factor)
    ] * (3 * n_segments)

    d_val = cuda.to_device(np.zeros(len(d_idx_i), dtype=np.float64))
    d_deri = cuda.to_device(np.zeros((n_segments, 3), dtype=np.float64))

    def cost_func(mu):
        return cuda_wrapper_chunked(
            mu,
            d_coords,
            d_times,
            d_idx_i,
            d_idx_j,
            sigma_nm,
            1.0,
            d_val,
            d_deri,
            chunk_size,
        )

    tmp = drift_est.reshape((-1, 3))
    for i in range(3):
        tmp[:, i] = convolve(tmp[:, i], np.ones(boxcar_width) / boxcar_width)
    drift_est = tmp.flatten()

    result = minimize(
        cost_func,
        drift_est,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"gtol": 1e-5, "ftol": 1e-12, "disp": False},
    )

    return result.x.reshape((n_segments, 3))


def temporal_refined_drift(
    dataset,
    max_drift_nm=300,
    segmentation_mode=1,
    segmentation_var=250,
    max_locs_per_segment=None,
    display=False,
):
    loc_frames = dataset[:, -1].astype(int)
    loc_coords = dataset[:, :3]

    base_result = segmentation_wrapper(
        loc_frames,
        segmentation_var,
        segmentation_mode,
        max_locs_per_segment,
        return_param_dict=True
    )
    base_segments = base_result.loc_segments
    base_center_frames = base_result.center_frames
    total_segments = len(base_center_frames)

    segmentations = []
    seg_count = total_segments
    while seg_count >= 3:
        segmentations.append(seg_count)
        seg_count = seg_count // 2
    segmentations = list(reversed(segmentations))

    valid_mask = base_result.loc_valid
    dataset[:, -1] = base_segments
    dataset = dataset[valid_mask]
    valid_segments = base_segments[valid_mask]
    loc_frames_valid = dataset[:, -1].astype(int)
    loc_coords_valid = dataset[:, :3].astype(np.float32)

    # Compute and push to GPU once
    idx_i, idx_j, successful = pair_indices_kdtree(loc_coords_valid, max_drift_nm)
    if not successful:
        raise MemoryError("pair_indices_kdtree ran out of memory; lower max_locs_per_segment "
                          "or max_drift_nm.")
    d_coords = cuda.to_device(loc_coords_valid)
    d_times = cuda.to_device(loc_frames_valid.astype(np.int32))
    d_idx_i = cuda.to_device(idx_i.astype(np.int32))
    d_idx_j = cuda.to_device(idx_j.astype(np.int32))

    common_frame_axis = np.arange(loc_frames_valid.min(), loc_frames_valid.max() + 1)

    # Assume initial sigma as one third of max drift
    sigma_nm = max_drift_nm / 3.0
    target_sigma_nm = 1.0

    last_valid_drift = None
    last_valid_center = None
    last_interpolated = None
    gradient_deltas = []

    if display:
        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    while sigma_nm >= target_sigma_nm:
        for i, n_segments in enumerate(segmentations):
            stride = total_segments // n_segments
            segment_ids = valid_segments // stride
            center_frames = [np.mean(loc_frames_valid[segment_ids == j]) for j in range(n_segments)]
            center_frames = np.array(center_frames)

            if last_valid_drift is not None:
                drift_init = interpolate_drift(last_valid_center, last_valid_drift, center_frames)
            else:
                drift_init = np.zeros((n_segments, 3))

            drift_est = optimize_drift(
                locs_nm=dataset,
                loc_frames=loc_frames_valid,
                n_segments=n_segments,
                sigma_nm=sigma_nm,
                drift_init=drift_init,
                drift_max_nm=max_drift_nm,
                d_coords=d_coords,
                d_times=d_times,
                d_idx_i=d_idx_i,
                d_idx_j=d_idx_j,
            )

            interpolated = interpolate_drift(center_frames, drift_est, common_frame_axis)
            interpolated -= np.mean(interpolated, axis=0)  # subtract mean for shift invariance

            if display:
                for dim, ax in enumerate(axes):
                    ax.plot(common_frame_axis, interpolated[:, dim], label=f"{n_segments} segments")

            if last_interpolated is not None:
                delta = compute_gradient_difference(interpolated, last_interpolated)
                gradient_deltas.append(delta)

                if i >= 2 and delta > gradient_deltas[-2]:
                    break

            last_valid_drift = drift_est
            last_valid_center = center_frames
            last_interpolated = interpolated

        if display:
            for dim, ax in enumerate(axes):
                ax.set_ylabel(["Drift X", "Drift Y", "Drift Z"][dim])
                ax.legend()
                ax.grid(True)
            axes[-1].set_xlabel("Frame")
            fig.suptitle(f"Drift Estimates at Sigma = {sigma_nm:.1f} nm", fontsize=14)
            plt.tight_layout()
            plt.show()

        sigma_nm /= 1.5  # refine sigma

    final_interp = interpolate_drift(last_valid_center, last_valid_drift, common_frame_axis)
    return final_interp, common_frame_axis


if __name__ == "__main__":
    h5_file = "../../data/mt_data_sim_2d_003.h5"
    dataset, gt_drift = load_simulation_dataset_and_gt_drift(h5_file)

    frames_per_segment = 10

    drift_cuda, _ = temporal_refined_drift(
        dataset=dataset,
        segmentation_mode=2,
        segmentation_var=frames_per_segment,
        max_locs_per_segment=50,
        max_drift_nm=300,
        display=True
    )
