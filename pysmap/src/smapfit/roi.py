"""Cutting fit ROIs around candidate positions."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .detect import Candidates


@dataclass
class ROIStack:
    """Square ROIs cut around candidates, ready for the fitter."""

    images: np.ndarray  # (n, roisize, roisize) float32, C-contiguous
    candidates: Candidates  # the candidate each ROI belongs to
    roisize: int
    mirrored: bool = False  # ROIs were flipped in x to match the PSF model

    def __len__(self) -> int:
        return self.images.shape[0]

    @property
    def half(self) -> int:
        return (self.roisize - 1) // 2

    def to_image_x(self, x_in_roi: np.ndarray) -> np.ndarray:
        """Map a fitted x inside the ROI back to image coordinates."""
        x = np.asarray(x_in_roi, dtype=np.float64)
        if self.mirrored:
            x = (self.roisize - 1) - x
        return x - self.half + self.candidates.x

    def to_image_y(self, y_in_roi: np.ndarray) -> np.ndarray:
        return np.asarray(y_in_roi, dtype=np.float64) - self.half + self.candidates.y


def cut_rois(photons: np.ndarray, candidates: Candidates, roisize: int,
             first_frame: int = 0, mirror: bool = False) -> ROIStack:
    """Cut ``roisize`` x ``roisize`` ROIs centred on each candidate.

    Candidates whose ROI would extend past the image border are dropped -- as
    ``simplefitter_cspline`` does.  SMAP's workflow instead shifts such ROIs
    inward and fits them anyway, which biases those localizations because the
    emitter is then no longer near the ROI centre.

    ``first_frame`` is the absolute frame number of ``photons[0]``, so that the
    absolute frame numbers carried by the candidates can be mapped onto this
    block.

    ``mirror`` flips the ROIs in x.  It is needed for calibrations that were
    computed from mirrored bead images (``parameters.emmirror`` in the
    ``_3Dcal.mat``); the flip is undone when the fitted position is mapped back
    to image coordinates, so the fitter itself never deals with mirroring.
    """
    if roisize % 2 == 0:
        raise ValueError(f"roisize must be odd, got {roisize}")

    photons = np.asarray(photons, dtype=np.float32)
    if photons.ndim == 2:
        photons = photons[None]
    n_frames, height, width = photons.shape
    half = (roisize - 1) // 2

    x, y = candidates.x, candidates.y
    inside = ((x >= half) & (x < width - half)
              & (y >= half) & (y < height - half))
    kept = candidates[inside]
    if len(kept) == 0:
        return ROIStack(np.empty((0, roisize, roisize), np.float32), kept,
                        roisize, mirror)

    span = np.arange(-half, half + 1)
    rows = kept.y[:, None] + span  # (n, roisize)
    cols = kept.x[:, None] + span
    frame_index = _block_index(kept.frame, first_frame, n_frames)

    images = photons[frame_index[:, None, None], rows[:, :, None], cols[:, None, :]]
    if mirror:
        images = images[:, :, ::-1]
    return ROIStack(np.ascontiguousarray(images, dtype=np.float32), kept,
                    roisize, mirror)


def _block_index(frame: np.ndarray, first_frame: int, n_frames: int) -> np.ndarray:
    """Absolute frame numbers -> position within the current block."""
    index = frame - first_frame
    if index.size and (index.min() < 0 or index.max() >= n_frames):
        raise ValueError(
            f"candidate frames {int(frame.min())}..{int(frame.max())} do not fit "
            f"the block starting at {first_frame} with {n_frames} frames")
    return index
