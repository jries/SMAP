"""Post-install smoke test: ``comet_self_test``.

Simulates a dataset with a known drift, runs COMET on it, and checks that the
recovered drift matches the ground truth. The backend is auto-detected so the
test works on a plain CPU-only machine as well as on a CUDA workstation.
"""

import argparse
import sys
import time
import traceback

import numpy as np

from comet.core.backends import best_backend, describe_backends
from comet.core.drift_optimizer import comet_run_kd

# Agreement with the ground-truth drift, averaged over the three axes.
THRESHOLD_NM = 1.0

# Pair count grows roughly quadratically with the localization count, so the
# NumPy backend needs a much smaller problem to finish in smoke-test time:
# ~3 s at the quick size versus several minutes at the full size.
QUICK_SIZE = dict(n_points=200, T=50, locs_per_frame=30)
FULL_SIZE = dict(n_points=1000, T=500, locs_per_frame=120)


def simulate_dataset(
    n_points=1000,         # base template molecules
    T=500,                 # timesteps / frames
    locs_per_frame=120,    # localizations per frame
    seed=123,
):
    rng = np.random.default_rng(seed)

    # Base coordinates (nm), spread in a 1000 nm cube
    gt_coords = rng.random((n_points, 3)) * 1000.0  # (N,3)

    # Time axis
    timesteps = np.arange(T)

    # 1D base drift waveform (nm / frame); then create 3D with cheap phase/scale offsets
    base = (3*np.sin(timesteps / 20.0 + 0.37) + 1.25 * np.cos(timesteps / 6.0 + 1.2) + timesteps/25) * 10.0
    dx = base
    dy = 0.7 * base
    dz = -0.5 * base

    # Per-frame drift (nm)
    gt_drift_nm = np.column_stack([dx, dy, dz])  # (T,3)

    # Vectorized localization generation
    M = T * locs_per_frame
    frame_idx = np.repeat(np.arange(T, dtype=np.int32), locs_per_frame)
    rng.shuffle(frame_idx)

    base_idx = rng.integers(0, n_points, size=M)
    base_xyz = gt_coords[base_idx]        # (M,3)
    drift_xyz = gt_drift_nm[frame_idx]    # (M,3)
    xyz = base_xyz + drift_xyz            # (M,3)

    # Pack into (N,4): x,y,z,frame
    locs = np.empty((M, 4), dtype=np.float32)
    locs[:, :3] = xyz.astype(np.float32)
    locs[:, 3] = frame_idx.astype(np.float32)
    return locs, gt_drift_nm


def run_mock_comet(dataset, display=False, mode=None):
    if mode is None:
        mode = best_backend()
    drift, _ = comet_run_kd(
        dataset=dataset,
        segmentation_mode=2,
        segmentation_var=2,
        initial_sigma_nm=120,
        max_drift_nm=100,
        target_sigma_nm=10,
        drift_max_bound_factor=2,
        boxcar_width=1,
        return_corrected_locs=True,
        interpolation_method='cubic',
        mode=mode,
        display=display
    )
    return drift[:, :3]


def drift_error_nm(estimated, ground_truth):
    """Offset-invariant agreement: mean per-axis std of the residual, in nm."""
    m = min(len(ground_truth), len(estimated))
    if m < 5:
        raise AssertionError("Too few frames for a meaningful check (need >= 5).")
    return float(np.mean(np.std(estimated[:m] - ground_truth[:m], axis=0)))


def _plot(gt, est):
    import matplotlib.pyplot as plt

    m = min(len(gt), len(est))
    fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    labels = ("x", "y", "z")
    for i, ax in enumerate(axes):
        ax.plot(gt[:m, i], label="GT {}".format(labels[i]))
        ax.plot(est[:m, i], "--", label="EST {}".format(labels[i]))
        ax.set_ylabel("nm")
        ax.legend()
    axes[-1].set_xlabel("Frame")
    fig.suptitle("COMET self-test: GT vs EST drift")
    plt.tight_layout()
    plt.show()


def main(argv=None):
    parser = argparse.ArgumentParser(description="COMET post-install self-test")
    parser.add_argument("--plot", action="store_true", help="Plot GT vs estimated drift")
    parser.add_argument("--mode", choices=["cuda", "torch", "cpu"], default=None,
                        help="Backend to test (default: auto-detect the fastest available)")
    parser.add_argument("--size", choices=["quick", "full"], default=None,
                        help="Problem size (default: full on a CUDA GPU, quick otherwise)")
    parser.add_argument("--verbose", action="store_true", help="Show optimizer progress")
    args = parser.parse_args(argv)

    import comet

    print("== COMET self-test ==")
    print("COMET  : {}".format(comet.__version__))
    print("Python : {}".format(sys.version.split()[0]))
    print("NumPy  : {}".format(np.__version__))
    for line in describe_backends():
        print("  " + line)

    mode = args.mode or best_backend()
    # Only a real CUDA GPU is fast enough to make the full size a sensible
    # default for a post-install check; MPS and CPU stay on the quick size.
    size = args.size or ("full" if mode == "cuda" else "quick")
    params = QUICK_SIZE if size == "quick" else FULL_SIZE
    print("Running on backend: {} ({} size, {} localizations)".format(
        mode, size, params["T"] * params["locs_per_frame"]))

    try:
        # Warm the backend on a tiny problem first, so JIT compilation
        # (numba caches the CPU kernel, but the CUDA kernel recompiles every
        # process) is not counted in the reported time.
        warm, _ = simulate_dataset(n_points=50, T=10, locs_per_frame=20)
        run_mock_comet(warm, display=False, mode=mode)

        dataset, gt = simulate_dataset(**params)
        t0 = time.perf_counter()
        est = run_mock_comet(dataset, display=args.verbose, mode=mode)
        elapsed = time.perf_counter() - t0
        err_nm = drift_error_nm(est, gt)
    except AssertionError as ae:
        print("FAILED - shape/consistency error: {}".format(ae))
        return 3
    except KeyboardInterrupt:
        print("Interrupted.")
        return 130
    except Exception as e:
        print("FAILED - unhandled error: {}".format(e))
        traceback.print_exc()
        return 1

    print("Drift correction took {:.2f} s on backend '{}' (excluding JIT warm-up)".format(
        elapsed, mode))
    print("Drift agreement with ground truth: {:.2f} nm".format(err_nm))

    if args.plot:
        _plot(gt, est)

    if np.isfinite(err_nm) and err_nm < THRESHOLD_NM:
        print("PASSED")
        return 0
    print("FAILED: {:.2f} nm >= {} nm threshold".format(err_nm, THRESHOLD_NM))
    return 2


if __name__ == "__main__":
    sys.exit(main())
