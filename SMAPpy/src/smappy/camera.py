"""Conversion of raw camera counts into photons."""

from __future__ import annotations

import numpy as np

from .metadata import CameraMetadata


def to_photons(frames: np.ndarray, cam: CameraMetadata) -> np.ndarray:
    """Convert raw camera counts to photons.

    ``frames`` may be a single image or a block of images ``(n, y, x)``.
    The result is float32.

    EM excess noise is *not* applied here.  The fitter assumes pure Poisson
    statistics; the EMCCD excess-noise factor is applied around the fit
    (see :attr:`CameraMetadata.excess_noise`) so that the fitter itself stays
    free of camera-specific assumptions.
    """
    cam.require("conversion", "offset")
    out = np.asarray(frames, dtype=np.float32) - np.float32(cam.offset)
    out *= np.float32(cam.adu_to_photons)
    return out
