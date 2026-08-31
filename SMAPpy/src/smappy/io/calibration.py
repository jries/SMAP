"""Reading of SMAP ``*_3Dcal.mat`` cubic-spline PSF calibration files.

This is the one module that has to understand MATLAB's array conventions.
Everything downstream works in a single, explicit convention:

* images are ``(y, x)`` C-contiguous, as they come out of a TIFF reader
* the spline coefficient array is stored as ``(64, nz, ny, nx)`` C-contiguous,
  so the flat index used by the fitter is ``((i*nz + z)*ny + y)*nx + x``
* fit parameters are ordered ``(x, y, photons, background, z)``

In the MATLAB fitter the equivalent array is ``(n1, n2, nz, 64)`` in
column-major order and the first axis -- called ``spline_xsize`` in the C
sources -- is the *first image axis*, i.e. rows, i.e. ``y``.  That is why
MATLAB's ``P(:,1)`` is y and ``P(:,2)`` is x.  We swap the two spatial axes
once, here, so that the rest of the code can call x "x".

The 64 coefficients per grid point are ordered ``i = 16*pz + 4*py + px``
(powers of the z, y, x fractional offsets), matching ``kernel_computeDelta3D``
in ``CPUsplineLib.cpp``; that ordering is preserved by the axis swap because
we swap the corresponding coefficient axes too.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.io


@dataclass
class SplineCalibration:
    """A cubic-spline PSF model loaded from a SMAP ``_3Dcal.mat`` file."""

    coeff: np.ndarray  # (64, nz, ny, nx) float32, C-contiguous
    dz: float  # z spacing of the calibration stack, nm
    z0: float  # index of the focal plane in the spline z grid
    x0: Optional[float] = None  # lateral centre of the calibration ROI, pixels
    psf: Optional[np.ndarray] = None  # (nz, ny, nx) averaged bead stack, if present
    em_on: Optional[bool] = None  # EM gain used for the bead stack
    em_mirror: Optional[bool] = None  # bead images were mirrored during calibration
    pixelsize_um: Optional[np.ndarray] = None
    source: Optional[Path] = None
    parameters: dict = field(default_factory=dict)  # raw calibration parameters

    @property
    def shape(self) -> tuple:
        """Spatial shape of the spline grid, ``(nz, ny, nx)``."""
        return self.coeff.shape[1:]

    def z_index_to_nm(self, z_index: np.ndarray) -> np.ndarray:
        """Convert a spline z index to nm relative to the focal plane."""
        return -(np.asarray(z_index) - self.z0) * self.dz

    def z_nm_to_index(self, z_nm: np.ndarray) -> np.ndarray:
        return -np.asarray(z_nm) / self.dz + self.z0


def _as_struct(obj):
    """Unwrap the 0-d object arrays scipy hands back for MATLAB structs."""
    while isinstance(obj, np.ndarray) and obj.dtype == object and obj.size == 1:
        obj = obj.item()
    return obj


def _first(obj):
    """First element of a MATLAB struct array (SXY may be 2-D for spatial cal)."""
    if isinstance(obj, np.ndarray) and obj.dtype == object and obj.size > 1:
        raise NotImplementedError(
            f"spatially varying calibration with {obj.shape} entries is not "
            "supported; only single-channel, single-region calibrations are"
        )
    return _as_struct(obj)


def _to_scalar(v, default=None):
    if v is None:
        return default
    v = np.asarray(v)
    return default if v.size == 0 else v.reshape(-1)[0]


def load_spline_calibration(path) -> SplineCalibration:
    """Load a SMAP ``_3Dcal.mat`` file (MATLAB v7 or v7.3)."""
    path = Path(path)
    mat = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)

    params = mat.get("parameters")
    if "SXY" in mat:
        entry = _first(mat["SXY"])
        cspline = _as_struct(entry.cspline)
        psf = getattr(entry, "PSF", None)
        em_on = _to_scalar(getattr(entry, "EMon", None))
    elif "cspline" in mat:
        # older flat layout written by calibrate3D_g for single-region fits
        entry = None
        cspline = _as_struct(mat["cspline"])
        psf = mat.get("PSF")
        em_on = None
    else:
        raise ValueError(f"{path.name}: no 'SXY' or 'cspline' found; not a 3D "
                         "calibration file")

    coeff = _as_struct(cspline.coeff)
    if isinstance(coeff, np.ndarray) and coeff.dtype == object:
        coeff = coeff.flat[0]  # cell array: single channel -> first cell
    coeff = np.asarray(coeff, dtype=np.float32)
    if coeff.ndim != 4 or coeff.shape[3] != 64:
        raise ValueError(f"{path.name}: unexpected coefficient shape "
                         f"{coeff.shape}, expected (n1, n2, nz, 64)")

    coeff = _matlab_coeff_to_c(coeff)

    if psf is not None:
        psf = np.asarray(psf, dtype=np.float32)
        if psf.ndim == 3:
            psf = np.ascontiguousarray(psf.transpose(2, 0, 1))  # -> (nz, ny, nx)
        else:
            psf = None

    par = {}
    if params is not None:
        params = _as_struct(params)
        for name in getattr(params, "_fieldnames", []):
            par[name] = getattr(params, name)

    em_mirror = _to_scalar(par.get("emmirror"))
    pixelsize = par.get("pixelsize")
    if pixelsize is not None:
        pixelsize = np.asarray(_as_struct(np.asarray(pixelsize).reshape(-1)[0]
                                          if np.asarray(pixelsize).dtype == object
                                          else pixelsize), dtype=float).reshape(-1)

    return SplineCalibration(
        coeff=coeff,
        dz=float(_to_scalar(cspline.dz)),
        z0=float(_to_scalar(cspline.z0)),
        x0=_maybe_float(_to_scalar(getattr(cspline, "x0", None))),
        psf=psf,
        em_on=None if em_on is None else bool(em_on),
        em_mirror=None if em_mirror is None else bool(em_mirror),
        pixelsize_um=pixelsize,
        source=path,
        parameters=par,
    )


def _maybe_float(v):
    return None if v is None else float(v)


def _matlab_coeff_to_c(coeff: np.ndarray) -> np.ndarray:
    """``(n1, n2, nz, 64)`` MATLAB order -> ``(64, nz, ny, nx)`` C order.

    MATLAB's first spatial axis is the image row axis (y), so the two spatial
    axes are swapped.  The 64 coefficients are indexed ``16*pz + 4*p1 + p2``
    where ``p1`` belongs to the first spatial axis; swapping the spatial axes
    therefore requires swapping the ``p1``/``p2`` digits of that index as well.
    """
    n1, n2, nz, _ = coeff.shape
    # (n1, n2, nz, 64) -> (64, nz, n1, n2), with n1 = y, n2 = x
    out = np.ascontiguousarray(coeff.transpose(3, 2, 0, 1))
    # reorder the coefficient axis: i = 16*pz + 4*py + px, with py from axis n1
    idx = np.arange(64)
    pz, p1, p2 = idx // 16, (idx // 4) % 4, idx % 4
    out = out[16 * pz + 4 * p2 + p1]
    return np.ascontiguousarray(out, dtype=np.float32)


def compute_delta3d(dx: float, dy: float, dz: float) -> np.ndarray:
    """The 64 monomials ``dz**pz * dy**py * dx**px``, ordered as the coefficients.

    Mirrors ``kernel_computeDelta3D``; ``px`` is the fastest-varying digit.
    """
    powers = np.arange(4)
    cz = np.asarray(dz, dtype=np.float64) ** powers
    cy = np.asarray(dy, dtype=np.float64) ** powers
    cx = np.asarray(dx, dtype=np.float64) ** powers
    return (cz[:, None, None] * cy[None, :, None] * cx[None, None, :]).reshape(64)


def evaluate_spline(cal: SplineCalibration, x: float, y: float, z: float,
                    roisize: int) -> np.ndarray:
    """Render the normalized PSF model into a ``(roisize, roisize)`` image.

    ``x``/``y`` are the emitter position in ROI pixel coordinates and ``z`` is
    in spline index units.  This reproduces the model evaluation of
    ``kernel_splineMLEFit_z_EMCCD`` and exists so the loaded coefficients can be
    checked against the bead stack stored in the same file.
    """
    _, nz, ny, nx = cal.coeff.shape

    xc = -((x - roisize / 2) + 0.5)
    yc = -((y - roisize / 2) + 0.5)
    xstart, ystart, zstart = int(np.floor(xc)), int(np.floor(yc)), int(np.floor(z))
    delta = compute_delta3d(xc - xstart, yc - ystart, z - zstart)

    # SMAP centres the (smaller) ROI inside the (larger) spline grid
    off_x = int(np.floor((nx + 1 - roisize) / 2))
    off_y = int(np.floor((ny + 1 - roisize) / 2))

    ix = np.clip(np.arange(roisize) + xstart + off_x, 0, nx - 1)
    iy = np.clip(np.arange(roisize) + ystart + off_y, 0, ny - 1)
    iz = int(np.clip(zstart, 0, nz - 1))

    block = cal.coeff[:, iz][:, iy[:, None], ix[None, :]]  # (64, roisize, roisize)
    return np.tensordot(delta.astype(np.float32), block, axes=(0, 0))


def warn_on_em_mismatch(cal: SplineCalibration, data_em_on: Optional[bool]) -> None:
    """Warn if bead calibration and data were taken with different EM settings."""
    if cal.em_on is None or data_em_on is None:
        return
    if bool(cal.em_on) != bool(data_em_on):
        warnings.warn(
            f"EM gain mismatch: calibration was acquired with EM "
            f"{'on' if cal.em_on else 'off'} but the data with EM "
            f"{'on' if data_em_on else 'off'}. The PSF model may be mirrored "
            "relative to the data. Acquire both with the same EM setting.",
            stacklevel=2,
        )
