"""Redundant cross-correlation drift estimation.

The other estimator in this package (:mod:`smappy.drift`) maximises the
overlap of individual localizations.  This one works on *images*: bin the
acquisition into time windows, render each one, and measure how far one has
moved against another by cross-correlation.  It is the classical method --
SMAP's ``finddriftfeature.m``, and Wang et al. 2014 for the redundancy -- and it
is here for two reasons: as an independent second opinion on a drift curve, and
as a coarse pass that costs seconds.

"Redundant" is the important part.  Correlating only consecutive windows would
accumulate every measurement error into the running sum.  Instead *every pair*
of windows is correlated, giving ``T(T-1)/2`` measurements of ``d_k - d_l`` for
``T`` unknowns, and the drift is the least-squares solution of that
overdetermined system.  A pair whose correlation peak was found in the wrong
place is then outvoted rather than believed, which is what the robust
reweighting below makes explicit.

Differences from SMAP's implementation, all deliberate:

* The correlation peak is located by a **quadratic fit** over a small window
  rather than a full 2D Gaussian fit.  The peak of a cross-correlation is not
  Gaussian anyway, both are accurate to a small fraction of a pixel, and this
  one is a linear least-squares problem with no starting values to get wrong.
* The overdetermined system is solved by **iteratively reweighted least
  squares** with a Cauchy weight, in place of MATLAB's ``nlinfit`` with robust
  options.  Same intent, and it is linear -- ``d_k - d_l`` is linear in the
  unknowns -- so it needs no nonlinear fit at all.
* Axial drift is **one-dimensional correlations of z histograms, computed in
  narrow slices of the field of view and summed** -- SMAP's
  ``finddisplacementZ``.  Correlating one pooled x-z image instead, which is
  what this module did first, mixes the z profiles of everything at different
  x and y before correlating: the local structure the peak comes from is gone,
  and the axial estimate is worthless (43 nm from the overlap estimator's, and
  barely better than no correction at all on the image).  Slice, correlate,
  then sum -- the sum of correlations is not the correlation of sums.
* Per-frame interpolation is **PCHIP**, which does not overshoot between time
  points the way a cubic spline through noisy estimates does.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .drift import Drift, _grouped, _positions_nm, _selection
from .locs import Localizations


@dataclass
class RCCSettings:
    """How finely to bin, render and search."""

    n_timepoints: int = 20           # time windows to correlate
    pixelsize_nm: float = 15.0       # rendering pixel size, lateral
    z_pixelsize_nm: float = 5.0      # z histogram bin width
    max_drift_nm: float = 500.0      # how far from zero shift the peak is sought
    fit_window: int = 3              # half-width of the sub-pixel fit, in pixels
    max_pixels: int = 2048           # cap on a rendered image, per axis
    use_z: Optional[bool] = None     # None: whenever the table has z_nm
    # Estimate from grouped localizations, as `DriftSettings` does by default,
    # so that the two estimators can be compared on the same footing.
    group: bool = True
    group_dx_nm: float = 50.0
    group_dt: int = 1
    # Axial: z histograms are correlated within tiles of the field of view and
    # the correlations summed, so that each one sees a *local* z profile.  The
    # tiles must be small in both directions -- with slices spanning y (SMAP's
    # `finddisplacementZ`, 200 nm in x only) each profile pools hundreds of
    # structures into a featureless slab and the correlation peak has nothing to
    # lock onto: measured on the clathrin data, the pairwise shifts came out
    # *uncorrelated* with the overlap estimator's (-0.09), against +0.96 with
    # 200 nm square tiles (SMAP's `finddisplacementZ2` tiles in y as well).
    tile_nm: float = 200.0
    tile_y_nm: Optional[float] = None    # None: square tiles, `tile_nm` in both
    # The axial peak is broad -- a z profile is far less structured than an
    # image -- so its search range is kept separate from, and smaller than, the
    # lateral one; a wide range lets the maximum wander hundreds of nm.
    axial_max_drift_nm: float = 200.0
    # Leave the zero-lag sample of the axial correlation out of the peak search
    # and the fit.  It carries anything that sits at the same z in both windows
    # regardless of drift -- a molecule on in both, or localizations pinned to a
    # particular z by the fit -- which pulls the estimate towards zero drift.
    exclude_zero_lag: bool = True


def estimate_drift_rcc(locs: Localizations, settings: Optional[RCCSettings] = None,
                       select=None, pixelsize_nm: Optional[float] = None,
                       display: bool = False) -> Drift:
    """Estimate drift by correlating every pair of time windows.

    Takes and returns the same things as :func:`smappy.drift.estimate_drift`,
    so the two are interchangeable and their results directly comparable.
    """
    settings = settings or RCCSettings()

    all_frames = np.asarray(locs["frame"], dtype=np.int64)
    n_frames = int(all_frames.max()) + 1 if len(all_frames) else 0

    source = locs[_selection(locs, select)]
    if settings.group:
        source = _grouped(source, settings, pixelsize_nm)
    x, y, z = _positions_nm(source, pixelsize_nm, settings.use_z)
    frames = np.asarray(source["frame"], dtype=np.int64)

    edges = np.linspace(0, n_frames, settings.n_timepoints + 1)
    window = np.clip(np.searchsorted(edges, frames, side="right") - 1,
                     0, settings.n_timepoints - 1)
    centres = (edges[:-1] + edges[1:]) / 2.0

    occupied = np.flatnonzero(np.bincount(window, minlength=settings.n_timepoints))
    if occupied.size < 2:
        raise ValueError("fewer than two time windows contain localizations")

    lateral = _pair_shifts(x, y, window, occupied, settings.pixelsize_nm,
                           settings, display, "x-y")
    dx = _solve(lateral[0], occupied.size)
    dy = _solve(lateral[1], occupied.size)

    if z is not None and np.any(z):
        # The axial pass slices the field of view, so it has to see the same
        # structure in each time window: take the lateral drift out first, or
        # a slice's contents differ between windows by the lateral drift and
        # the z profiles being correlated are not of the same thing.
        lateral_per_frame = _to_frames(centres[occupied],
                                       np.column_stack([dx, dy, np.zeros_like(dx)]),
                                       n_frames)
        dz = _solve(_axial_shifts(x - lateral_per_frame[frames, 0],
                                  y - lateral_per_frame[frames, 1],
                                  z, window, occupied, settings, display),
                    occupied.size)
    else:
        dz = np.zeros(occupied.size)

    per_window = np.column_stack([dx, dy, dz])
    drift = _to_frames(centres[occupied], per_window, n_frames)
    return Drift(drift, None, n_used=len(x))


# ------------------------------------------------------------------- correlating
def _render(a, b, window, occupied, grid) -> np.ndarray:
    """One 2D histogram per time window, all on the same grid."""
    (a0, a1, na), (b0, b1, nb) = grid
    images = np.zeros((occupied.size, nb, na), dtype=np.float32)
    for k, w in enumerate(occupied):
        inside = window == w
        images[k] = np.histogram2d(b[inside], a[inside], bins=(nb, na),
                                   range=((b0, b1), (a0, a1)))[0]
    return images


def _grid(values, pixelsize: float, max_pixels: int):
    """Axis extent and bin count, folded back if the field of view is too wide.

    Folding (SMAP does the same) keeps the FFT bounded on a large field: the
    structure repeats, but a shift of the whole image is still a shift.
    """
    lo, hi = float(np.min(values)), float(np.max(values))
    n = int(round((hi - lo) / pixelsize)) + 1
    if n > max_pixels:
        span = max_pixels * pixelsize
        values = lo + np.mod(values - lo, span)
        lo, hi, n = float(np.min(values)), float(np.max(values)), max_pixels
    return values, (lo, hi, n)


def _pair_shifts(a, b, window, occupied, pixelsize, settings: RCCSettings,
                 display: bool, label: str):
    """Shift of every window against every other, in nm, as two (T, T) arrays."""
    a, grid_a = _grid(a, pixelsize, settings.max_pixels)
    b, grid_b = _grid(b, pixelsize, settings.max_pixels)
    images = _render(a, b, window, occupied, (grid_a, grid_b))

    # The correlation is circular, so the padding has to leave room for every
    # lag that will be looked at: without it a short axis (an x-z image is only
    # a few tens of pixels tall) wraps its own signal into the search window and
    # the shift comes out meaningless.
    reach = max(2, int(round(settings.max_drift_nm / pixelsize)))
    shape = tuple(1 << int(np.ceil(np.log2(max(s + reach + 1, 64))))
                  for s in images.shape[1:])
    spectra = np.fft.rfft2(images, s=shape)

    n = occupied.size
    shifts = np.zeros((2, n, n))
    for k in range(n):
        for l in range(k + 1, n):
            correlation = np.fft.fftshift(
                np.fft.irfft2(spectra[k] * np.conj(spectra[l]), s=shape))
            row, col = _peak(correlation, reach, settings.fit_window)
            shifts[0, k, l], shifts[0, l, k] = col * pixelsize, -col * pixelsize
            shifts[1, k, l], shifts[1, l, k] = row * pixelsize, -row * pixelsize
    if display:
        print(f"RCC: {n} time windows, {n * (n - 1) // 2} {label} correlations "
              f"on {shape[1]}x{shape[0]} pixels")
    return shifts


def _axial_shifts(x, y, z, window, occupied, settings: RCCSettings, display: bool):
    """Shift in z of every time window against every other, in nm.

    One z histogram per (time window, slice of the field of view); the
    correlation of a pair is summed over the slices, so every slice contributes
    its own local z profile instead of everything being pooled first.
    """
    z_pixel = settings.z_pixelsize_nm
    z_lo, z_hi = float(np.min(z)), float(np.max(z))
    n_bins = max(8, int(round((z_hi - z_lo) / z_pixel)) + 1)
    z_index = np.clip(((z - z_lo) / z_pixel).astype(np.int64), 0, n_bins - 1)

    tile = ((x - x.min()) / settings.tile_nm).astype(np.int64)
    column = ((y - y.min()) / (settings.tile_y_nm or settings.tile_nm)).astype(np.int64)
    tile = tile * (column.max() + 1) + column
    _, tile = np.unique(tile, return_inverse=True)
    n_tiles = int(tile.max()) + 1

    position = np.full(int(occupied.max()) + 1, -1, dtype=np.int64)
    position[occupied] = np.arange(occupied.size)
    slot = position[window]

    n = occupied.size
    flat = (slot * n_tiles + tile) * n_bins + z_index
    histograms = np.bincount(flat, minlength=n * n_tiles * n_bins).astype(np.float32)
    histograms = histograms.reshape(n, n_tiles, n_bins)
    # subtracting each histogram's mean removes the pedestal that would
    # otherwise dominate the correlation over its structure
    histograms -= histograms.mean(axis=2, keepdims=True)

    reach = max(2, int(round(settings.axial_max_drift_nm / z_pixel)))
    size = 1 << int(np.ceil(np.log2(n_bins + reach + 1)))
    spectra = np.fft.rfft(histograms, n=size, axis=2).astype(np.complex64)
    if display:
        print(f"RCC: {n} time windows, {n * (n - 1) // 2} axial correlations over "
              f"{n_tiles} tiles of {settings.tile_nm:.0f} nm, {n_bins} z bins")

    shifts = np.zeros((n, n))
    for k in range(n):
        for l in range(k + 1, n):
            correlation = np.fft.fftshift(
                np.fft.irfft((spectra[k] * np.conj(spectra[l])).sum(axis=0), n=size))
            lag = _peak_1d(correlation, reach, settings.fit_window,
                           settings.exclude_zero_lag)
            shifts[k, l], shifts[l, k] = lag * z_pixel, -lag * z_pixel
    return shifts


def _peak_1d(correlation: np.ndarray, reach: int, half: int,
             exclude_zero: bool = True) -> float:
    """Sub-pixel lag of a 1D correlation maximum, ignoring the zero lag.

    Zero lag carries a spike from molecules that are on in both windows -- the
    same emitter at the same z -- which says nothing about drift, so it is left
    out of both the search and the fit, as SMAP does.
    """
    centre = correlation.size // 2
    lo, hi = max(0, centre - reach), min(correlation.size, centre + reach + 1)
    patch = correlation[lo:hi].copy()
    if exclude_zero:
        patch[centre - lo] = -np.inf
    index = int(np.argmax(patch)) + lo

    left, right = index - half, index + half + 1
    if left < 0 or right > correlation.size:
        return float(index - centre)
    lags = np.arange(left, right) - centre
    values = correlation[left:right].astype(float)
    keep = (lags != 0) if exclude_zero else np.ones(lags.size, bool)
    if keep.sum() < 3:
        return float(index - centre)
    a, b, _ = np.polyfit(lags[keep], values[keep], 2)
    if a >= 0:                        # not a maximum: keep the sampled lag
        return float(index - centre)
    offset = -b / (2 * a)
    return offset if abs(offset - (index - centre)) <= half else float(index - centre)


def _peak(correlation: np.ndarray, reach: int, half: int) -> Tuple[float, float]:
    """Sub-pixel position of the correlation maximum, relative to zero shift."""
    ny, nx = correlation.shape
    cy, cx = ny // 2, nx // 2
    y0, y1 = max(0, cy - reach), min(ny, cy + reach + 1)
    x0, x1 = max(0, cx - reach), min(nx, cx + reach + 1)
    patch = correlation[y0:y1, x0:x1]

    # smooth before picking the pixel, as SMAP does: the maximum of a noisy
    # correlation can sit one pixel off the ridge, and the fit below is local
    smoothed = _boxcar(patch, 5)
    row, col = np.unravel_index(int(np.argmax(smoothed)), patch.shape)
    row, col = row + y0, col + x0

    if not (half <= row < ny - half and half <= col < nx - half):
        return float(row - cy), float(col - cx)
    offset = _quadratic_peak(correlation[row - half:row + half + 1,
                                         col - half:col + half + 1])
    return float(row - cy) + offset[0], float(col - cx) + offset[1]


def _boxcar(image: np.ndarray, width: int) -> np.ndarray:
    from scipy.ndimage import uniform_filter

    return uniform_filter(image, width, mode="nearest")


def _quadratic_peak(patch: np.ndarray) -> Tuple[float, float]:
    """Offset of the maximum of a 2D quadratic fitted to ``patch``, in pixels."""
    half = patch.shape[0] // 2
    v, u = np.mgrid[-half:half + 1, -half:half + 1]
    u, v = u.ravel().astype(float), v.ravel().astype(float)
    design = np.column_stack([np.ones_like(u), u, v, u * u, u * v, v * v])
    c, *_ = np.linalg.lstsq(design, patch.ravel().astype(float), rcond=None)

    # solve the 2x2 system for the stationary point of the quadratic
    hessian = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])
    try:
        du, dv = np.linalg.solve(hessian, [-c[1], -c[2]])
    except np.linalg.LinAlgError:
        return 0.0, 0.0
    if not (abs(du) <= half and abs(dv) <= half):  # fit ran away: keep the pixel
        return 0.0, 0.0
    return float(dv), float(du)


# ----------------------------------------------------------------------- solving
def _solve(shifts: np.ndarray, n: int, iterations: int = 5) -> np.ndarray:
    """Least-squares drift per window from all pairwise shifts, robustly.

    Every pair contributes one equation ``d_k - d_l = shift``; the system is
    overdetermined by construction, which is the whole point of measuring every
    pair.  Reweighting then lets a pair whose correlation peak was misplaced be
    outvoted instead of believed.  The drift is fixed to zero mean, since only
    differences are measured.
    """
    k, l = np.triu_indices(n, 1)
    rows = np.arange(k.size)
    design = np.zeros((k.size + 1, n))
    design[rows, k] = 1.0
    design[rows, l] = -1.0
    design[-1, :] = 1.0                      # gauge: the drift averages to zero
    observed = np.append(shifts[k, l], 0.0)

    weights = np.ones(observed.size)
    solution = np.zeros(n)
    for _ in range(iterations):
        scaled = design * weights[:, None]
        solution, *_ = np.linalg.lstsq(scaled, observed * weights, rcond=None)
        residual = design @ solution - observed
        scale = 1.4826 * np.median(np.abs(residual[:-1] - np.median(residual[:-1])))
        if scale <= 0:
            break
        weights = 1.0 / (1.0 + (residual / (3.0 * scale)) ** 2)  # Cauchy
        weights[-1] = 1.0                    # never downweight the gauge
    return solution


def _to_frames(centres: np.ndarray, per_window: np.ndarray, n_frames: int) -> np.ndarray:
    """Per-window drift to per-frame, without overshooting between windows."""
    from scipy.interpolate import PchipInterpolator

    frames = np.arange(n_frames)
    if len(centres) < 2:
        return np.repeat(per_window[:1], n_frames, axis=0)
    return np.column_stack([
        PchipInterpolator(centres, per_window[:, axis], extrapolate=True)(frames)
        for axis in range(3)])
