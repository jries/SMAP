"""Check a loaded spline calibration against the bead stack in the same file.

The check is non-circular: the spline coefficients and the averaged
experimental PSF stack are independent arrays inside the ``_3Dcal.mat``.  With
zero fractional offset the cubic spline evaluates to its zeroth coefficient, so
``coeff[0]`` must reproduce the bead stack at the corresponding knots.  The
transposed / mirrored / z-flipped comparisons are controls: if the axis
convention in ``calibration.py`` were wrong, one of those would win instead.

Usage:  python scripts/check_calibration.py <file_3dcal.mat>
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from smappy.io.calibration import load_spline_calibration, evaluate_spline  # noqa: E402


def main(path: str) -> int:
    cal = load_spline_calibration(path)
    print(f"coeff  {cal.coeff.shape}  dz={cal.dz} nm  z0={cal.z0}  x0={cal.x0}")
    print(f"PSF    {None if cal.psf is None else cal.psf.shape}"
          f"  EMon={cal.em_on}  emmirror={cal.em_mirror}")
    if cal.psf is None:
        print("no PSF stack stored in this file: cannot check the axis convention")
        return 1

    c0 = cal.coeff[0]
    nz, ny, nx = c0.shape
    # the spline grid is shorter than the stack by 3 points per axis; the knot
    # of grid point i is stack sample i+1
    ref = cal.psf[1:1 + nz, 1:1 + ny, 1:1 + nx]
    good = np.isfinite(ref)

    # the spline is normalized differently from the (max-normalized) stack
    scale = np.sum(c0[good] * ref[good]) / np.sum(ref[good] ** 2)
    resid = np.abs(c0[good] - scale * ref[good])

    def corr(a, b):
        g = np.isfinite(b)
        return np.corrcoef(a[g], b[g])[0, 1]

    print(f"\nnormalization coeff/PSF: {scale:.5f}")
    print(f"residual after scaling:  max={resid.max() / np.nanmax(ref):.5f} "
          f"(relative)  rms={np.sqrt((resid ** 2).mean()):.6f}")
    print(f"\ncorrelation, as loaded : {corr(c0, ref):.6f}   <- must be the largest")
    print(f"  control, transposed  : {corr(c0, np.ascontiguousarray(ref.transpose(0, 2, 1))):.6f}")
    print(f"  control, x-mirrored  : {corr(c0, np.ascontiguousarray(ref[:, :, ::-1])):.6f}")
    print(f"  control, z-flipped   : {corr(c0, np.ascontiguousarray(ref[::-1])):.6f}")

    # the fitter's own parametrization must agree with the raw coefficients
    r, z = 13, int(round(cal.z0))
    model = evaluate_spline(cal, r / 2 - 0.5, r / 2 - 0.5, float(z), r)
    off_x = int(np.floor((nx + 1 - r) / 2))
    off_y = int(np.floor((ny + 1 - r) / 2))
    dev = np.abs(model - c0[z, off_y:off_y + r, off_x:off_x + r]).max()
    print(f"\nevaluate_spline at a knot vs coeff[0]: max deviation {dev:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
