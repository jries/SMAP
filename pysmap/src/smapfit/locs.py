"""Localizations: fit results as a table.

Positions are kept in **camera pixels**, which is what the fitter measures.
Converting to nm needs the pixel size, a separate calibration that can change
or be corrected afterwards, so it is a derived quantity produced by
:func:`to_nm` at the end of the pipeline rather than baked into the fit.

Two things are undone here because they are properties of the camera, not
units: the ROI offset on the chip is added to the pixel coordinates, so
positions from different acquisition ROIs are comparable, and photon counts are
scaled back by the EMCCD excess-noise factor the fitter was shielded from.

z is the exception: it comes from the calibration's ``dz`` in nm and has no
pixel equivalent, so ``z_nm`` is always in nm.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterator, Optional

import numpy as np

from .metadata import CameraMetadata
from .psf import FitResult, PSFModel
from .roi import ROIStack


@dataclass
class Localizations:
    """A table of localizations: named columns of equal length."""

    columns: Dict[str, np.ndarray] = field(default_factory=dict)
    metadata: Dict[str, object] = field(default_factory=dict)

    def __len__(self) -> int:
        return 0 if not self.columns else len(next(iter(self.columns.values())))

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.columns[key]
        return Localizations({k: v[key] for k, v in self.columns.items()},
                             dict(self.metadata))

    def __contains__(self, key: str) -> bool:
        return key in self.columns

    def __iter__(self) -> Iterator[str]:
        return iter(self.columns)

    def keys(self):
        return self.columns.keys()

    def extend(self, other: "Localizations") -> None:
        """Append another table with the same columns."""
        if not self.columns:
            self.columns = {k: np.asarray(v) for k, v in other.columns.items()}
            return
        if set(self.columns) != set(other.columns):
            raise ValueError(
                f"column mismatch: {sorted(set(self.columns) ^ set(other.columns))}")
        for k in self.columns:
            self.columns[k] = np.concatenate([self.columns[k], other.columns[k]])

    def __str__(self) -> str:
        return (f"{len(self)} localizations, columns: "
                f"{', '.join(sorted(self.columns))}")


def fit_to_localizations(result: FitResult, rois: ROIStack, model: PSFModel,
                         cam: CameraMetadata) -> Localizations:
    """Convert raw fit output into a localization table, in camera pixels."""
    p = model.unpack(result)
    excess = cam.excess_noise
    roi_x, roi_y = cam.roi_offset

    # ROI pixel -> image pixel (this undoes the mirror, if one was applied)
    # -> chip pixel, by adding the offset of the acquisition ROI
    x_pix = rois.to_image_x(p["x_roi"]) + roi_x
    y_pix = rois.to_image_y(p["y_roi"]) + roi_y

    cols: Dict[str, np.ndarray] = {
        "frame": rois.candidates.frame.astype(np.int64),
        "x_pix": x_pix,
        "y_pix": y_pix,
        "photons": p["photons"] * excess,
        "background": p["background"] * excess,
        "x_err_pix": p["x_err_pix"],
        "y_err_pix": p["y_err_pix"],
        "photons_err": p["photons_err"] * excess,
        "background_err": p["background_err"] * excess,
        "logl": result.logl,
        # log-likelihood per pixel, comparable between ROI sizes and cameras
        "logl_rel": result.logl * excess / rois.roisize ** 2,
        "peak_x_pix": rois.candidates.x + roi_x,
        "peak_y_pix": rois.candidates.y + roi_y,
        "iterations": result.iterations.astype(np.int32),
    }

    if "z_nm" in p:  # z has no pixel equivalent: it comes from the calibration
        cols["z_nm"] = p["z_nm"]
        cols["z_err_nm"] = p["z_err_nm"]
    for name in ("sigma_pix", "sigma_x_pix", "sigma_y_pix"):
        if name in p:
            cols[name] = p[name]

    cols["loc_precision_pix"] = np.sqrt((cols["x_err_pix"] ** 2
                                         + cols["y_err_pix"] ** 2) / 2)

    for key, value in cols.items():
        if key not in ("frame", "iterations"):
            cols[key] = np.asarray(value, dtype=np.float32)

    return Localizations(cols, {"units": "pixel"})


# columns in pixels, and the nm name they get; z is already in nm
_PIXEL_COLUMNS = {
    "x_pix": "x_nm", "y_pix": "y_nm",
    "x_err_pix": "x_err_nm", "y_err_pix": "y_err_nm",
    "peak_x_pix": "peak_x_nm", "peak_y_pix": "peak_y_nm",
    "loc_precision_pix": "loc_precision_nm",
    "sigma_pix": "sigma_nm", "sigma_x_pix": "sigma_x_nm",
    "sigma_y_pix": "sigma_y_nm",
}


def to_nm(locs: Localizations, pixelsize_nm: float,
          keep_pixels: bool = False) -> Localizations:
    """Return a copy with pixel columns converted to nm.

    ``keep_pixels`` keeps the pixel columns alongside the nm ones.  ``z_nm`` is
    already in nm and is left untouched.
    """
    if not locs.columns:
        return locs
    if locs.metadata.get("units") == "nm":
        return locs

    columns = {}
    for name, values in locs.columns.items():
        target = _PIXEL_COLUMNS.get(name)
        if target is None:
            columns[name] = values
            continue
        columns[target] = np.asarray(values, dtype=np.float32) * np.float32(pixelsize_nm)
        if keep_pixels:
            columns[name] = values

    metadata = dict(locs.metadata)
    metadata["units"] = "pixel+nm" if keep_pixels else "nm"
    metadata["pixelsize_nm"] = float(pixelsize_nm)
    return Localizations(columns, metadata)


def valid(locs: Localizations, max_fit_distance: Optional[float] = None
          ) -> np.ndarray:
    """Mask of usable localizations.

    Drops non-finite results, and -- if ``max_fit_distance`` is given (in
    pixels) -- fits that ran away from their candidate position.  Since the
    fitter no longer clamps x and y to the ROI centre, this is where such fits
    are rejected rather than silently parked on a boundary.
    """
    mask = np.ones(len(locs), dtype=bool)
    for key in ("x_pix", "y_pix", "photons"):
        mask &= np.isfinite(locs[key])
    if "z_nm" in locs:
        mask &= np.isfinite(locs["z_nm"])
    if max_fit_distance is not None:
        mask &= (np.abs(locs["x_pix"] - locs["peak_x_pix"]) < max_fit_distance)
        mask &= (np.abs(locs["y_pix"] - locs["peak_y_pix"]) < max_fit_distance)
    return mask
