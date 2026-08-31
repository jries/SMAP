"""Axis-convention tests for the spline calibration reader.

These need a real SMAP _3Dcal.mat; set SMAPPY_TEST_CAL to point at one.
"""
import os
import numpy as np
import pytest

from smappy.io.calibration import (
    compute_delta3d, evaluate_spline, load_spline_calibration,
)

CAL = os.environ.get("SMAPPY_TEST_CAL")
needs_cal = pytest.mark.skipif(not CAL, reason="set SMAPPY_TEST_CAL")


def test_compute_delta3d_ordering():
    d = compute_delta3d(2.0, 3.0, 5.0)
    for pz in range(4):
        for py in range(4):
            for px in range(4):
                assert d[16 * pz + 4 * py + px] == pytest.approx(5.0**pz * 3.0**py * 2.0**px)
    assert compute_delta3d(0.0, 0.0, 0.0)[0] == 1.0
    assert compute_delta3d(0.0, 0.0, 0.0)[1:].sum() == 0.0


@needs_cal
def test_coefficients_match_bead_stack():
    """coeff[0] is the model value at a knot, so it must match the PSF stack."""
    cal = load_spline_calibration(CAL)
    assert cal.psf is not None
    c0 = cal.coeff[0]
    nz, ny, nx = c0.shape
    ref = cal.psf[1:1 + nz, 1:1 + ny, 1:1 + nx]
    g = np.isfinite(ref)

    def corr(b):
        gb = np.isfinite(b)
        return np.corrcoef(c0[gb], b[gb])[0, 1]

    direct = corr(ref)
    assert direct > 0.99
    # wrong axis conventions must be clearly worse
    assert direct > corr(np.ascontiguousarray(ref.transpose(0, 2, 1))) + 0.1
    assert direct > corr(np.ascontiguousarray(ref[:, :, ::-1])) + 0.1
    assert direct > corr(np.ascontiguousarray(ref[::-1])) + 0.1


@needs_cal
def test_evaluate_spline_reproduces_coefficients_at_knot():
    cal = load_spline_calibration(CAL)
    _, nz, ny, nx = cal.coeff.shape
    r, z = 13, int(round(cal.z0))
    model = evaluate_spline(cal, r / 2 - 0.5, r / 2 - 0.5, float(z), r)
    ox = int(np.floor((nx + 1 - r) / 2))
    oy = int(np.floor((ny + 1 - r) / 2))
    assert np.allclose(model, cal.coeff[0, z, oy:oy + r, ox:ox + r])
