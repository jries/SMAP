"""PSF models: what to fit, and how to read the fitted parameters.

A model knows how to call the fitter and how to turn the raw parameter vector
into named quantities in physical units.  Everything else in the pipeline is
model-agnostic, so adding a model means adding a class here (and a kernel in
``csrc/models.hpp``), not touching the pipeline.

Fitted parameters are always ordered ``(x, y, photons, background, ...)`` with
x the ROI column and y the row.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

from . import _fit3d
from .io.calibration import SplineCalibration


@dataclass
class FitResult:
    """Raw fitter output for a stack of ROIs."""

    theta: np.ndarray  # (n, NV) fitted parameters
    crlb: np.ndarray  # (n, NV) Cramer-Rao lower bounds (variances)
    logl: np.ndarray  # (n,) log-likelihood
    iterations: np.ndarray  # (n,) iterations used

    def __len__(self) -> int:
        return self.theta.shape[0]


class PSFModel:
    """Base class; see :class:`GaussianPSF` and :class:`SplinePSF`."""

    name = "psf"

    def fit(self, rois: np.ndarray, iterations: int = 50,
            n_threads: int = 0) -> FitResult:
        raise NotImplementedError

    def unpack(self, result: FitResult) -> Dict[str, np.ndarray]:
        """Named quantities derived from the raw parameters."""
        raise NotImplementedError

    @property
    def is_3d(self) -> bool:
        return False


@dataclass
class GaussianPSF(PSFModel):
    """Gaussian PSF with a free width.

    ``elliptical=False`` fits one width (``sigma``); ``elliptical=True`` fits
    ``sigma_x`` and ``sigma_y`` separately, which is what an astigmatic
    calibration measurement needs.
    """

    sigma: float = 1.0
    elliptical: bool = False
    name = "gauss"

    def fit(self, rois: np.ndarray, iterations: int = 50,
            n_threads: int = 0) -> FitResult:
        fn = _fit3d.fit_gauss_xy if self.elliptical else _fit3d.fit_gauss_free
        return FitResult(*fn(rois, float(self.sigma), int(iterations),
                             int(n_threads)))

    def unpack(self, result: FitResult) -> Dict[str, np.ndarray]:
        p, crlb = result.theta, np.clip(result.crlb, 0, None)
        out = {
            "x_roi": p[:, 0], "y_roi": p[:, 1],
            "photons": p[:, 2], "background": p[:, 3],
            "x_err_pix": np.sqrt(crlb[:, 0]), "y_err_pix": np.sqrt(crlb[:, 1]),
            "photons_err": np.sqrt(crlb[:, 2]),
            "background_err": np.sqrt(crlb[:, 3]),
        }
        if self.elliptical:
            out["sigma_x_pix"] = p[:, 4]
            out["sigma_y_pix"] = p[:, 5]
        else:
            out["sigma_pix"] = p[:, 4]
            out["sigma_err_pix"] = np.sqrt(crlb[:, 4])
        return out

    def __str__(self) -> str:
        kind = "elliptical" if self.elliptical else "free sigma"
        return f"GaussianPSF({kind}, start sigma={self.sigma:g} pix)"


@dataclass
class SplinePSF(PSFModel):
    """Experimental PSF interpolated by a cubic spline, from a SMAP calibration.

    ``z_start`` is the starting z in nm relative to focus; it is converted to
    the spline's index units.  ``mirror`` reports whether the calibration was
    computed from mirrored bead images, in which case the ROIs must be flipped
    in x before fitting (see :func:`smappy.roi.cut_rois`).
    """

    calibration: SplineCalibration
    z_start_nm: float = 0.0
    name = "cspline"

    @property
    def is_3d(self) -> bool:
        return True

    @property
    def mirror(self) -> bool:
        """Whether ROIs must be x-flipped to match this calibration."""
        return bool(self.calibration.em_mirror)

    def fit(self, rois: np.ndarray, iterations: int = 50,
            n_threads: int = 0) -> FitResult:
        z_start = float(self.calibration.z_nm_to_index(self.z_start_nm))
        return FitResult(*_fit3d.fit_cspline(rois, self.calibration.coeff,
                                             z_start, int(iterations),
                                             int(n_threads)))

    def unpack(self, result: FitResult) -> Dict[str, np.ndarray]:
        p, crlb = result.theta, np.clip(result.crlb, 0, None)
        cal = self.calibration
        return {
            "x_roi": p[:, 0], "y_roi": p[:, 1],
            "photons": p[:, 2], "background": p[:, 3],
            "z_nm": cal.z_index_to_nm(p[:, 4]),
            "x_err_pix": np.sqrt(crlb[:, 0]), "y_err_pix": np.sqrt(crlb[:, 1]),
            "photons_err": np.sqrt(crlb[:, 2]),
            "background_err": np.sqrt(crlb[:, 3]),
            "z_err_nm": np.sqrt(crlb[:, 4]) * cal.dz,
        }

    def __str__(self) -> str:
        nz, ny, nx = self.calibration.shape
        return (f"SplinePSF({nx}x{ny}x{nz} grid, dz={self.calibration.dz:g} nm"
                f"{', mirrored calibration' if self.mirror else ''})")
