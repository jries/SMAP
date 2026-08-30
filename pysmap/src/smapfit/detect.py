"""Candidate detection: image filtering and local-maximum finding.

The image is filtered to suppress noise and (for the difference-of-Gaussians)
the smooth background, then every strict local maximum above a cutoff becomes a
candidate.  Both stages work on a block of frames ``(n, y, x)`` at once.

Differences to SMAP's ``ImageFilter`` / ``PeakFinder``:

* no Anscombe transform (it belongs with a background estimate, which is not
  part of this pipeline yet) and no sCMOS variance weighting
* no background offset subtracted before filtering.  In SMAP that value only
  compensated the zero padding of the convolution at the image borders; we pad
  by edge replication instead, which removes the artifact at its source
* only the difference-of-Gaussians and Gaussian filters, and only the dynamic
  and absolute cutoffs
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Protocol

import numpy as np
from scipy import ndimage

try:
    from . import _fit3d
except ImportError:  # pure-Python fallback if the extension is not built
    _fit3d = None


# --------------------------------------------------------------------- filters
class ImageFilter(Protocol):
    """Turns a block of images into a block of detection images."""

    def __call__(self, images: np.ndarray, n_threads: int = 0) -> np.ndarray: ...


@dataclass(frozen=True)
class GaussFilter:
    """Gaussian smoothing; use when the background is flat or already removed."""

    sigma: float = 1.2
    use_extension: bool = True

    def __call__(self, images: np.ndarray, n_threads: int = 0) -> np.ndarray:
        if self.sigma <= 0:
            return np.asarray(images, dtype=np.float32)
        radius = _gauss_radius(self.sigma)
        if self.use_extension and _fit3d is not None:
            return _fit3d.filter_images(
                _as_block(images), self.sigma, 0.0, radius,
                separable=radius >= _SEPARABLE_FROM_RADIUS, n_threads=n_threads)
        return _gaussian(images, self.sigma, radius)

    def __str__(self) -> str:
        return f"Gauss(sigma={self.sigma:g})"


@dataclass(frozen=True)
class DoGFilter:
    """Difference of Gaussians: smooths the PSF and removes smooth background.

    The wide Gaussian follows SMAP: ``max(1, 2.5 * sigma)``.  The kernel sums to
    zero exactly, so a constant background cancels.
    """

    sigma: float = 1.2
    use_extension: bool = True

    @property
    def sigma_wide(self) -> float:
        return max(1.0, 2.5 * self.sigma)

    def __call__(self, images: np.ndarray, n_threads: int = 0) -> np.ndarray:
        radius = _dog_radius(self.sigma)
        if self.use_extension and _fit3d is not None:
            # small kernel: one 2-D convolution with the subtracted kernels
            # (the 2-D DoG has rank 2, so it is not separable, but it can be
            # applied directly).  Large kernel: two Gaussians fused into one
            # separable traversal.
            return _fit3d.filter_images(
                _as_block(images), self.sigma, self.sigma_wide, radius,
                separable=radius >= _SEPARABLE_FROM_RADIUS, n_threads=n_threads)
        narrow = _gaussian(images, self.sigma, radius)
        wide = _gaussian(images, self.sigma_wide, radius)
        return narrow - wide

    def __str__(self) -> str:
        return f"DoG(sigma={self.sigma:g}, wide={self.sigma_wide:g})"


def _as_block(images: np.ndarray) -> np.ndarray:
    """A C-contiguous float32 block ``(n, y, x)``, as the extension expects."""
    out = np.ascontiguousarray(images, dtype=np.float32)
    return out[None] if out.ndim == 2 else out


def _gauss_radius(sigma: float) -> int:
    """Kernel radius, following SMAP's ``max(5, ceil(3.5/2*s)*2+1)`` window."""
    return max(2, int(np.ceil(1.75 * sigma)))


def _dog_radius(sigma: float) -> int:
    """Kernel radius, following SMAP's ``max(ceil(6*s-1), 3)`` window."""
    window = max(int(np.ceil(6.0 * sigma - 1.0)), 3)
    return max(1, (window - 1) // 2)


# Below this radius a single 2-D convolution beats the separable passes: it
# does more arithmetic per pixel ((2r+1)^2 vs 4(2r+1)) but in one traversal
# with no temporaries, and for small kernels the traversal is what costs.
# Measured crossover on an M1 Pro; see scripts/benchmark.py.
_SEPARABLE_FROM_RADIUS = 4


def _gaussian(images: np.ndarray, sigma: float, radius: int) -> np.ndarray:
    """Separable Gaussian over the two image axes, edge-replicate padding."""
    images = np.asarray(images, dtype=np.float32)
    axes = (images.ndim - 2, images.ndim - 1)
    out = images
    for axis in axes:
        out = ndimage.gaussian_filter1d(out, sigma, axis=axis, mode="nearest",
                                        radius=radius, output=np.float32)
    return out


# --------------------------------------------------------------------- cutoffs
@dataclass(frozen=True)
class AbsoluteCutoff:
    """Fixed threshold on the filtered image."""

    value: float

    def __call__(self, maxima_values: np.ndarray) -> float:
        return float(self.value)

    def __str__(self) -> str:
        return f"absolute({self.value:g})"


@dataclass(frozen=True)
class DynamicCutoff:
    """Threshold derived from the distribution of local-maximum intensities.

    Follows SMAP: the cutoff is the median of the maxima plus ``factor`` times
    the slope of their 20-80 % quantile range.  Computed per frame, so it
    adapts to changing background without needing the whole movie.
    """

    factor: float = 1.7

    def __call__(self, maxima_values: np.ndarray) -> float:
        n = maxima_values.size
        if n == 0:
            return 0.0
        if n < 10:
            return float(np.mean(maxima_values) * self.factor)
        q20, q50, q80 = np.quantile(maxima_values, (0.2, 0.5, 0.8))
        slope = (q80 - q20) / 0.6
        return float(q50 + slope * self.factor)

    def __str__(self) -> str:
        return f"dynamic(factor={self.factor:g})"


# ------------------------------------------------------------------ candidates
@dataclass
class Candidates:
    """Peak candidates found in a block of frames (pixel coordinates)."""

    frame: np.ndarray  # frame index, int64
    x: np.ndarray  # column, int32
    y: np.ndarray  # row, int32
    value: np.ndarray  # intensity in the filtered image, float32

    def __len__(self) -> int:
        return self.frame.size

    def __getitem__(self, mask) -> "Candidates":
        return Candidates(self.frame[mask], self.x[mask], self.y[mask],
                          self.value[mask])

    @classmethod
    def empty(cls) -> "Candidates":
        return cls(np.empty(0, np.int64), np.empty(0, np.int32),
                   np.empty(0, np.int32), np.empty(0, np.float32))


_HOLLOW3 = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=bool)


def local_maxima(images: np.ndarray, threshold: float = -np.inf,
                 use_extension: bool = True, n_threads: int = 0) -> tuple:
    """Strict 3x3 local maxima -> ``(frame, y, x, value)``.

    "Strict" means larger than all eight neighbours, as in SMAP's
    ``maximumfind2.c``; plateaus therefore yield no maximum.  The one-pixel
    border is excluded.  ``threshold`` skips pixels at or below it, which is
    only a speed-up: pass ``-inf`` (the default) when the cutoff still has to be
    derived from the distribution of all maxima.

    Uses the C++ single-pass search when the extension is available -- it beats
    a neighbourhood filter because short-circuit evaluation rejects almost every
    pixel after one or two comparisons.  The numpy path is kept as a reference
    implementation and produces identical results.
    """
    images = np.ascontiguousarray(images, dtype=np.float32)
    if images.ndim == 2:
        images = images[None]

    if use_extension and _fit3d is not None:
        frame, y, x, value = _fit3d.find_maxima(images, float(threshold),
                                                n_threads)
        return frame.astype(np.int64), y.astype(np.int64), x.astype(np.int64), value

    centre = images[:, 1:-1, 1:-1]
    is_max = ((centre > images[:, :-2, :-2]) & (centre > images[:, :-2, 1:-1])
              & (centre > images[:, :-2, 2:]) & (centre > images[:, 1:-1, :-2])
              & (centre > images[:, 1:-1, 2:]) & (centre > images[:, 2:, :-2])
              & (centre > images[:, 2:, 1:-1]) & (centre > images[:, 2:, 2:]))
    if np.isfinite(threshold):
        is_max &= centre > threshold
    frame, y, x = np.nonzero(is_max)
    y += 1
    x += 1
    return frame, y, x, images[frame, y, x]


def find_candidates(filtered: np.ndarray, cutoff,
                    n_threads: int = 0) -> Candidates:
    """All local maxima of ``filtered`` above the cutoff.

    The cutoff is evaluated per frame, from that frame's maxima, which is what
    makes :class:`DynamicCutoff` adaptive.
    """
    filtered = np.asarray(filtered, dtype=np.float32)
    if filtered.ndim == 2:
        filtered = filtered[None]

    # an absolute cutoff can be applied during the search, which is cheaper
    fixed = cutoff(np.empty(0)) if isinstance(cutoff, AbsoluteCutoff) else -np.inf
    frames, rows, cols, values = local_maxima(filtered, threshold=fixed,
                                              n_threads=n_threads)
    if frames.size == 0:
        return Candidates.empty()
    if np.isfinite(fixed):
        return Candidates(frames, cols.astype(np.int32), rows.astype(np.int32),
                          values)

    # maxima come out sorted by frame, so each frame is a contiguous slice
    keep = np.zeros(frames.size, dtype=bool)
    n_frames = filtered.shape[0]
    bounds = np.searchsorted(frames, np.arange(n_frames + 1))
    for start, stop in zip(bounds[:-1], bounds[1:]):
        if stop > start:
            block = values[start:stop]
            keep[start:stop] = block > cutoff(block)

    return Candidates(frame=frames[keep].astype(np.int64),
                      x=cols[keep].astype(np.int32),
                      y=rows[keep].astype(np.int32),
                      value=values[keep])


@dataclass
class PeakFinder:
    """Filter a block of frames and return the candidates above the cutoff."""

    filter: ImageFilter
    cutoff: object
    n_threads: int = 0  # 0 = one per core

    def __call__(self, photons: np.ndarray, first_frame: int = 0,
                 n_threads: Optional[int] = None):
        """Returns ``(candidates, filtered)``; frame indices are absolute."""
        threads = self.n_threads if n_threads is None else n_threads
        filtered = self.filter(photons, n_threads=threads)
        candidates = find_candidates(filtered, self.cutoff, n_threads=threads)
        candidates.frame = candidates.frame + first_frame
        return candidates, filtered

    def __str__(self) -> str:
        return f"PeakFinder({self.filter}, {self.cutoff})"
