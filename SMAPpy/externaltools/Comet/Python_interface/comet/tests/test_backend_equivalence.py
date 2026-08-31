"""Cross-backend agreement.

The cuda and torch backends reimplement the same cost function as the NumPy
reference, so they must agree with it. These tests skip automatically when the
backend is not available, which makes the file a no-op on a CPU-only machine.
"""

import numpy as np
import pytest

from comet.core.drift_optimizer import comet_run_kd
from comet.tests.conftest import drift_residual_nm

# Backends solve the same problem in different float orders, so allow a little
# more slack than the CPU-vs-ground-truth tolerance.
BACKEND_AGREEMENT_NM = 1.0

RUN_KWARGS = dict(
    segmentation_mode=2,
    segmentation_var=2,
    initial_sigma_nm=120,
    max_drift_nm=100,
    target_sigma_nm=10,
    boxcar_width=1,
    interpolation_method="cubic",
    display=False,
)


def run(dataset, mode):
    return comet_run_kd(dataset=dataset.copy(), mode=mode, **RUN_KWARGS)


@pytest.fixture(scope="module")
def cpu_reference(request):
    """CPU drift estimate for the shared dataset, computed once per module."""
    from comet.tests.conftest import make_drifting_dataset
    locs, gt_drift = make_drifting_dataset()
    return locs, gt_drift, run(locs, "cpu")


@pytest.mark.cuda
def test_cuda_matches_cpu(cpu_reference):
    locs, _, reference = cpu_reference
    assert drift_residual_nm(run(locs, "cuda"), reference) < BACKEND_AGREEMENT_NM


@pytest.mark.cuda
def test_cuda_matches_ground_truth(cpu_reference):
    locs, gt_drift, _ = cpu_reference
    assert drift_residual_nm(run(locs, "cuda"), gt_drift) < BACKEND_AGREEMENT_NM


@pytest.mark.torch
def test_torch_matches_cpu(cpu_reference):
    locs, _, reference = cpu_reference
    assert drift_residual_nm(run(locs, "torch"), reference) < BACKEND_AGREEMENT_NM


@pytest.mark.torch
def test_torch_matches_ground_truth(cpu_reference):
    locs, gt_drift, _ = cpu_reference
    assert drift_residual_nm(run(locs, "torch"), gt_drift) < BACKEND_AGREEMENT_NM


def test_every_available_backend_recovers_the_drift(any_backend):
    """Runs once per backend present on this machine."""
    from comet.tests.conftest import make_drifting_dataset
    locs, gt_drift = make_drifting_dataset()

    drift = run(locs, any_backend)

    assert np.isfinite(drift).all()
    assert drift_residual_nm(drift, gt_drift) < BACKEND_AGREEMENT_NM


@pytest.mark.slow
@pytest.mark.cuda
def test_cuda_on_the_full_test_dataset():
    """Regression check against the large checked-in dataset, if present."""
    import pathlib

    from comet.core.io_utils import load_thunderstorm_csv

    csv_file = pathlib.Path(__file__).resolve().parents[3] / "test_dataset" / "test_dataset.csv"
    if not csv_file.exists():
        pytest.skip("test_dataset/test_dataset.csv is not present in this checkout")

    dataset = load_thunderstorm_csv(str(csv_file))
    drift = comet_run_kd(
        dataset=dataset, segmentation_mode=2, segmentation_var=50,
        initial_sigma_nm=40, max_drift_nm=100, target_sigma_nm=10,
        drift_max_bound_factor=2, boxcar_width=1,
        interpolation_method="cubic", mode="cuda", display=False,
    )

    assert drift.shape[1] == 4
    assert np.isfinite(drift).all(), "drift output contains NaNs"
