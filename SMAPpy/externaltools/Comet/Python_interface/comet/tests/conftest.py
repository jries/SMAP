"""Shared pytest configuration and fixtures for the COMET test suite.

The suite is designed to pass on a plain CPU-only machine: anything needing a
GPU carries a marker and is skipped automatically when the backend is absent.
"""

import os

# Must happen before anything imports matplotlib.pyplot, otherwise a GUI
# backend may be selected and tests that plot will block or fail headless.
os.environ.setdefault("MPLBACKEND", "Agg")

import subprocess
import sys

import numpy as np
import pytest


def run_python(args, timeout=600):
    """Run the current interpreter with `args` and capture its output.

    Spelled the pre-3.7 way (capture_output= and text= were added in 3.7) so
    the suite also runs on the Python 3.6/3.7 legacy lane.
    """
    return subprocess.run(
        [sys.executable] + list(args),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, timeout=timeout,
    )

from comet.core.backends import cuda_available, torch_available

# Markers and the predicate that decides whether the backend is usable here.
_BACKEND_MARKERS = {
    "cuda": (cuda_available, "no NVIDIA GPU / CUDA driver available"),
    "torch": (torch_available, "torch is not installed"),
}


def pytest_configure(config):
    config.addinivalue_line("markers", "cuda: needs an NVIDIA GPU and a working CUDA driver")
    config.addinivalue_line("markers", "torch: needs the optional torch dependency")
    config.addinivalue_line("markers", "slow: takes more than a few seconds")


def pytest_addoption(parser):
    parser.addoption("--run-slow", action="store_true", default=False,
                     help="also run tests marked slow")


def pytest_runtest_setup(item):
    """Skip GPU-dependent and slow tests unless the environment supports them."""
    for marker in item.iter_markers():
        probe = _BACKEND_MARKERS.get(marker.name)
        if probe is not None:
            is_available, reason = probe
            if not is_available():
                pytest.skip(reason)
        if marker.name == "slow" and not item.config.getoption("--run-slow"):
            pytest.skip("slow test; pass --run-slow to enable")


def make_drifting_dataset(n_points=200, n_frames=50, locs_per_frame=30, seed=0,
                          drift_amplitude_nm=30.0, extent_nm=1000.0, jitter_nm=0.0):
    """Synthetic localizations with a known, smooth drift.

    A fixed set of template molecules is sampled every frame and displaced by a
    smooth per-frame drift, which is what COMET is meant to recover.

    Returns
    -------
    locs : ndarray of shape (n_frames * locs_per_frame, 4)
        Columns [x_nm, y_nm, z_nm, frame].
    gt_drift : ndarray of shape (n_frames, 3)
        The per-frame drift that was applied, in nm.
    """
    rng = np.random.default_rng(seed)
    template = rng.random((n_points, 3)) * extent_nm

    t = np.arange(n_frames)
    base = drift_amplitude_nm * (np.sin(t / 12.0) + t / n_frames)
    gt_drift = np.column_stack([base, 0.7 * base, -0.5 * base])

    frame_idx = np.repeat(t, locs_per_frame)
    pick = rng.integers(0, n_points, size=frame_idx.size)
    xyz = template[pick] + gt_drift[frame_idx]
    if jitter_nm:
        xyz = xyz + rng.normal(0.0, jitter_nm, xyz.shape)

    locs = np.empty((frame_idx.size, 4), dtype=np.float64)
    locs[:, :3] = xyz
    locs[:, 3] = frame_idx
    return locs, gt_drift


def drift_residual_nm(estimated, ground_truth):
    """Offset-invariant agreement between two drift traces, in nm.

    Drift is only defined up to a constant offset, so compare the standard
    deviation of the residual rather than its absolute value.
    """
    m = min(len(estimated), len(ground_truth))
    return float(np.mean(np.std(estimated[:m, :3] - ground_truth[:m, :3], axis=0)))


@pytest.fixture
def drifting_dataset():
    """Small synthetic dataset: (locs, gt_drift). Fast enough for the CPU backend."""
    return make_drifting_dataset()


@pytest.fixture
def thunderstorm_csv(tmp_path, drifting_dataset):
    """Write the synthetic dataset out as a 3D ThunderSTORM CSV and return the path."""
    locs, _ = drifting_dataset
    path = tmp_path / "synthetic_locs.csv"
    with open(str(path), "w") as fh:
        fh.write('"id","frame","x [nm]","y [nm]","z [nm]"\n')
        for i, row in enumerate(locs):
            fh.write("{},{},{:.3f},{:.3f},{:.3f}\n".format(
                i + 1, int(row[3]), row[0], row[1], row[2]))
    return path


@pytest.fixture(params=["cpu", "cuda", "torch"])
def any_backend(request):
    """Parametrize over every backend, skipping those unavailable here."""
    backend = request.param
    probe = _BACKEND_MARKERS.get(backend)
    if probe is not None and not probe[0]():
        pytest.skip(probe[1])
    return backend
