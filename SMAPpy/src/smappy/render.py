"""Rendering localizations into a superresolution image.

Two stages, deliberately separate:

* **accumulation** produces float planes -- a weight plane (localizations, or
  photons, per superresolution pixel) and, when localizations carry individual
  colours, three colour planes.  This is the expensive part and depends on the
  data and the field of view.
* **display** turns those planes into RGB: normalise, gamma, look-up table.
  This is cheap and depends only on the settings, so a viewer can change the
  contrast or the colour map without re-rendering.

SMAP splits the same way (``renderSMAP`` / ``drawerSMAP``).

Colour: an intensity image is recoloured at display time, one lookup per output
pixel.  Colouring by a field (z, frame, ...) cannot work that way, because two
localizations of different colour in one pixel have to mix, so each one
contributes its own colour during accumulation.  We keep **four** planes then --
the three colour-weighted ones plus the plain weight -- which allows both
composites at display time from a single render:

``"hue"``    ``colour_sum / weight_sum`` gives the hue, ``weight_sum`` the
             brightness.  Saturation and gamma then act on brightness only and
             leave the hue alone.
``"sum"``    SMAP's composite: display ``colour_sum`` directly.  Identical to
             ``"hue"`` wherever the pixel is below saturation and gamma is 1;
             above it, per-channel clipping pulls dense regions towards white.

Coordinates are in whatever unit the localizations are in (nm, normally).  A
pixel of the rendered image covers ``[x0 + j*p, x0 + (j+1)*p)``, so a position
maps to pixel ``floor((x - x0) / p)``.  SMAP instead rounds and shifts the range
by half a pixel beforehand; the result is the same grid, expressed once instead
of twice.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from . import lut as luts
from .locs import Localizations

try:
    from . import _render
except ImportError:  # pure-Python fallback if the extension is not built
    _render = None


@dataclass(frozen=True)
class FieldOfView:
    """The rendered grid: origin, pixel size and shape.

    ``x0``/``y0`` are the coordinates of the *outer edge* of pixel (0, 0).
    Pixels are square, which keeps the image undistorted under any zoom.
    """

    x0: float
    y0: float
    pixelsize: float
    nx: int
    ny: int

    @property
    def x1(self) -> float:
        return self.x0 + self.nx * self.pixelsize

    @property
    def y1(self) -> float:
        return self.y0 + self.ny * self.pixelsize

    @property
    def shape(self) -> Tuple[int, int]:
        return (self.ny, self.nx)

    @property
    def extent(self) -> Tuple[float, float, float, float]:
        """``(left, right, bottom, top)``, ready for ``imshow(extent=...)``."""
        return (self.x0, self.x1, self.y1, self.y0)

    @classmethod
    def from_range(cls, xrange, yrange, pixelsize: float) -> "FieldOfView":
        """Cover the given coordinate ranges at a fixed pixel size."""
        nx = max(int(round((xrange[1] - xrange[0]) / pixelsize)), 1)
        ny = max(int(round((yrange[1] - yrange[0]) / pixelsize)), 1)
        return cls(float(xrange[0]), float(yrange[0]), float(pixelsize), nx, ny)

    @classmethod
    def fit(cls, xrange, yrange, nx: int, ny: int) -> "FieldOfView":
        """Fit the ranges into a canvas of ``nx`` x ``ny`` pixels.

        The pixel size follows from the canvas, which is what a viewer wants:
        the cost of a render then depends on the window, not on the zoom.  The
        ranges are centred, and the larger one decides the pixel size so
        nothing is cut off.
        """
        pixelsize = max((xrange[1] - xrange[0]) / nx, (yrange[1] - yrange[0]) / ny)
        cx = 0.5 * (xrange[0] + xrange[1])
        cy = 0.5 * (yrange[0] + yrange[1])
        return cls(cx - 0.5 * nx * pixelsize, cy - 0.5 * ny * pixelsize,
                   float(pixelsize), int(nx), int(ny))

    @classmethod
    def around(cls, x, y, pixelsize: float, margin: float = 0.0) -> "FieldOfView":
        """The smallest field of view containing all of ``x``, ``y``."""
        return cls.from_range((np.min(x) - margin, np.max(x) + margin),
                              (np.min(y) - margin, np.max(y) + margin), pixelsize)

    def to_pixels(self, x, y) -> Tuple[np.ndarray, np.ndarray]:
        """Data coordinates -> fractional pixel coordinates."""
        inv = 1.0 / self.pixelsize
        return ((np.asarray(x, dtype=np.float64) - self.x0) * inv,
                (np.asarray(y, dtype=np.float64) - self.y0) * inv)

    def zoomed(self, factor: float, cx: Optional[float] = None,
               cy: Optional[float] = None) -> "FieldOfView":
        """Same canvas, pixel size scaled by ``factor`` about ``(cx, cy)``.

        ``factor < 1`` zooms in.  The centre defaults to the centre of the
        current view; passing the cursor position is what makes scroll-zoom
        feel right.
        """
        cx = 0.5 * (self.x0 + self.x1) if cx is None else cx
        cy = 0.5 * (self.y0 + self.y1) if cy is None else cy
        p = self.pixelsize * factor
        return FieldOfView(cx - (cx - self.x0) * factor,
                           cy - (cy - self.y0) * factor, p, self.nx, self.ny)

    def moved(self, dx: float, dy: float) -> "FieldOfView":
        return FieldOfView(self.x0 + dx, self.y0 + dy, self.pixelsize,
                           self.nx, self.ny)


@dataclass
class RenderedImage:
    """Accumulated planes plus the geometry they were rendered on."""

    fov: FieldOfView
    weight: np.ndarray                      # (ny, nx)
    color: Optional[np.ndarray] = None      # (ny, nx, 3), colour-weighted
    n_locs: int = 0

    @property
    def is_colored(self) -> bool:
        return self.color is not None


def render_histogram(x, y, fov: FieldOfView, weights=None, colors=None
                     ) -> RenderedImage:
    """Count localizations (or their weights) per pixel.

    ``colors`` is an ``(n, 3)`` array of per-localization RGB; pass it to get a
    field-coloured image.
    """
    px, py = fov.to_pixels(x, y)
    col = np.floor(px).astype(np.int64)
    row = np.floor(py).astype(np.int64)
    inside = (col >= 0) & (col < fov.nx) & (row >= 0) & (row < fov.ny)

    index = (row[inside] * fov.nx + col[inside])
    w = None if weights is None else np.asarray(weights, dtype=np.float64)[inside]
    n_pixels = fov.nx * fov.ny

    weight = np.bincount(index, weights=w, minlength=n_pixels)
    weight = weight.reshape(fov.shape).astype(np.float32)

    color = None
    if colors is not None:
        c = np.asarray(colors, dtype=np.float64)[inside]
        if w is not None:
            c = c * w[:, None]
        color = np.stack(
            [np.bincount(index, weights=c[:, k], minlength=n_pixels)
             for k in range(3)], axis=-1).reshape(fov.ny, fov.nx, 3).astype(np.float32)

    return RenderedImage(fov, weight, color, int(inside.sum()))


# SMAP's dynamic contrast, and its parameterisation: a single number p that
# saturates 10^-p of the pixels.  The value distribution of a superresolution
# image is extremely heavy-tailed -- a fiducial bead outweighs the structure by
# three orders of magnitude -- so a fixed maximum is useless and a quantile is
# the only scale that transfers between datasets.  p is the natural handle
# because the useful range spans decades: 2 is bright, 4 is dark.
DEFAULT_CONTRAST = 3.0


def saturated_fraction(contrast: float) -> float:
    """The fraction of pixels driven to full scale: 10^-p."""
    return 10.0 ** -float(contrast)


def normalize(image: np.ndarray, imax: Optional[float] = None,
              contrast: float = DEFAULT_CONTRAST) -> Tuple[np.ndarray, float]:
    """Scale to [0, 1] and clip.  Returns the image and the value used as 1.

    ``imax`` overrides the dynamic contrast with an absolute value.  As in SMAP,
    a quantile that lands on zero -- a very sparse image -- falls back to the
    maximum.
    """
    if imax is None:
        imax = float(np.quantile(image, 1.0 - saturated_fraction(contrast)))
        if imax <= 0:
            imax = float(image.max())
    if imax <= 0:
        return np.zeros_like(image), 0.0
    return np.clip(image / imax, 0.0, 1.0), imax


def to_rgb(rendered: RenderedImage, lut: luts.LUT = "hot", invert: bool = False,
           imax: Optional[float] = None, contrast: float = DEFAULT_CONTRAST,
           gamma: float = 1.0, color_mode: str = "hue") -> np.ndarray:
    """Turn accumulated planes into a float32 ``(ny, nx, 3)`` RGB image in [0, 1].

    For an intensity image ``lut`` recolours the normalised weight plane.  For a
    field-coloured one the LUT was already applied per localization, and
    ``color_mode`` selects the composite (see the module docstring).
    """
    if rendered.color is None:
        norm, _ = normalize(rendered.weight, imax, contrast)
        if gamma != 1.0:
            norm = norm ** gamma
        table = luts.get(lut, invert)
        index = np.clip((norm * len(table)).astype(np.int32), 0, len(table) - 1)
        return table[index]

    if color_mode == "sum":
        rgb, _ = normalize(rendered.color, imax, contrast)
        return rgb ** gamma if gamma != 1.0 else rgb
    if color_mode != "hue":
        raise ValueError(f"unknown color_mode {color_mode!r}; use 'hue' or 'sum'")

    brightness, _ = normalize(rendered.weight, imax, contrast)
    if gamma != 1.0:
        brightness = brightness ** gamma
    with np.errstate(invalid="ignore", divide="ignore"):
        hue = rendered.color / rendered.weight[..., None]
    hue[~np.isfinite(hue)] = 0.0
    return (hue * brightness[..., None]).astype(np.float32)


# --------------------------------------------------------------- sigma policy
@dataclass(frozen=True)
class SigmaSettings:
    """How a localization precision becomes a rendering sigma.

    Straight from SMAP, because it is what makes precision-weighted images
    readable rather than a field of single bright pixels: blow the precision up
    a little (``factor``), never let a blob fall below a floor -- an absolute
    one in data units and one tied to the rendered pixel size -- and cap the
    outliers, since a handful of localizations with an absurd precision would
    otherwise dominate the render time.
    """

    factor: float = 1.0
    min_sigma: float = 0.0          # data units, SMAP's mingaussnm
    min_sigma_pixels: float = 0.7   # rendered pixels, SMAP's mingausspix
    max_factor: float = 10.0        # cap at this multiple of the median

    def apply(self, precision, pixelsize: float) -> np.ndarray:
        """Localization precision -> rendering sigma, in data units."""
        sigma = np.asarray(precision, dtype=np.float32) * np.float32(self.factor)
        finite = sigma[np.isfinite(sigma)]
        cap = (self.max_factor * np.median(finite)) if finite.size else np.inf
        floor = max(self.min_sigma, self.min_sigma_pixels * pixelsize)
        return np.clip(np.nan_to_num(sigma, nan=floor), floor, cap).astype(np.float32)


# SMAP's roiks: the kernel is cut off at this many sigma.
ROI_SIGMA = 2.7


def render(x, y, fov: FieldOfView, sigma=None, weights=None, colors=None,
           roi_sigma: float = ROI_SIGMA, n_threads: int = 0,
           use_extension: bool = True) -> RenderedImage:
    """Accumulate localizations into ``fov``.

    ``sigma`` selects the mode and is in data units: ``None`` renders a
    histogram, a scalar renders every localization with the same Gaussian, and
    an array renders each with its own -- normally a localization precision put
    through :class:`SigmaSettings`.  There is one kernel behind all three; see
    ``csrc/render.hpp``.

    ``weights`` is the intensity per localization (one value broadcasts), and
    ``colors`` an ``(n, 3)`` array for field colouring.
    """
    x = np.ascontiguousarray(x, dtype=np.float32)
    y = np.ascontiguousarray(y, dtype=np.float32)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be matching 1-D arrays")

    gaussian = sigma is not None
    sigma_x = sigma_y = np.zeros(1, np.float32)
    if gaussian:
        sigma = np.ascontiguousarray(sigma, dtype=np.float32).reshape(-1)
        if sigma.size not in (1, x.size):
            raise ValueError("sigma must be a scalar or one value per localization")
        sigma_x = sigma_y = sigma

    w = np.ones(1, np.float32) if weights is None else \
        np.ascontiguousarray(weights, dtype=np.float32).reshape(-1)
    c = np.zeros((0, 3), np.float32) if colors is None else \
        np.ascontiguousarray(colors, dtype=np.float32).reshape(-1, 3)
    if c.size and c.shape[0] != x.size:
        raise ValueError("colors must have one row per localization")

    if use_extension and _render is not None:
        weight, color, n = _render.render(
            x, y, sigma_x, sigma_y, w, c, fov.x0, fov.y0, fov.pixelsize,
            fov.nx, fov.ny, gaussian, roi_sigma, 512, n_threads)
        return RenderedImage(fov, weight, color, int(n))

    return _render_numpy(x, y, fov, sigma_x if gaussian else None, w, c, roi_sigma)


def _render_numpy(x, y, fov, sigma, weights, colors, roi_sigma) -> RenderedImage:
    """Reference implementation of :func:`render`; the extension must match it."""
    from scipy.special import erf

    if sigma is None:
        w = None if weights.size == 1 else weights
        img = render_histogram(x, y, fov, weights=w,
                               colors=colors if colors.size else None)
        if weights.size == 1:
            img.weight *= weights[0]
            if img.color is not None:
                img.color *= weights[0]
        return img

    weight = np.zeros(fov.shape, np.float64)
    color = np.zeros((fov.ny, fov.nx, 3), np.float64) if colors.size else None
    px, py = fov.to_pixels(x, y)
    n_locs = 0

    def kernel(center, position, s):
        """Pixel-integrated Gaussian over the ROI, and its sum."""
        d = int(roi_sigma * s + 1.0)
        edges = (np.arange(-d, d + 2) - (position - center)) / (np.sqrt(2) * s)
        k = 0.5 * np.diff(erf(edges))
        return np.arange(center - d, center + d + 1), k, k.sum()

    for i in range(x.size):
        s = (sigma[0] if sigma.size == 1 else sigma[i]) / fov.pixelsize
        cx, cy = int(np.floor(px[i])), int(np.floor(py[i]))
        if 0 <= cx < fov.nx and 0 <= cy < fov.ny:
            n_locs += 1
        cols, kx, sx = kernel(cx, px[i], s)
        rows, ky, sy = kernel(cy, py[i], s)
        tile = np.outer(ky, kx) * (weights[0] if weights.size == 1 else weights[i]) \
            / (sx * sy)
        keep_r = (rows >= 0) & (rows < fov.ny)
        keep_c = (cols >= 0) & (cols < fov.nx)
        if not keep_r.any() or not keep_c.any():
            continue
        block = np.ix_(rows[keep_r], cols[keep_c])
        tile = tile[np.ix_(keep_r, keep_c)]
        weight[block] += tile
        if color is not None:
            color[block] += tile[..., None] * colors[i]

    return RenderedImage(fov, weight.astype(np.float32),
                         None if color is None else color.astype(np.float32), n_locs)


# ------------------------------------------------------------------ front door
# The columns a table can carry positions and precisions in, best first.  Both
# unit systems work; the field of view just has to be in the same one.
POSITION_FIELDS = (("x_nm", "y_nm"), ("x_pix", "y_pix"))
PRECISION_FIELDS = ("loc_precision_nm", "loc_precision_pix")


def _pick(locs: Localizations, names, what: str):
    for name in names:
        if (all(n in locs for n in name) if isinstance(name, tuple)
                else name in locs):
            return name
    raise KeyError(f"no {what} column in the table; looked for {names}")


def positions(locs: Localizations, select=None):
    """The position columns of a table, in whatever unit it carries."""
    x_name, y_name = _pick(locs, POSITION_FIELDS, "position")
    x, y = locs[x_name], locs[y_name]
    return (x, y) if select is None else (x[select], y[select])


@dataclass
class RenderSettings:
    """What to accumulate: the kernel, and what colours the localizations."""

    mode: str = "precision"     # "hist", "gauss" (one sigma), "precision"
    sigma: float = 10.0         # data units, for mode="gauss"
    sigma_settings: "SigmaSettings" = field(default_factory=lambda: SigmaSettings())
    roi_sigma: float = ROI_SIGMA
    color_field: Optional[str] = None    # None: an intensity image, recoloured later
    color_range: Optional[Tuple[float, float]] = None
    weight_field: Optional[str] = None   # None: one count per localization


@dataclass
class DisplaySettings:
    """How to turn accumulated planes into RGB.  Changing these re-displays, not
    re-renders -- except ``lut``/``invert`` when colouring by a field, because
    then the colour is baked into the accumulation."""

    lut: luts.LUT = "hot"
    invert: bool = False
    contrast: float = DEFAULT_CONTRAST   # saturate 10^-contrast of the pixels
    imax: Optional[float] = None         # an absolute scale, overriding contrast
    gamma: float = 1.0
    color_mode: str = "hue"

    def apply(self, rendered: RenderedImage) -> np.ndarray:
        return to_rgb(rendered, self.lut, self.invert, self.imax, self.contrast,
                      self.gamma, self.color_mode)


def render_locs(locs: Localizations, fov: FieldOfView,
                settings: Optional[RenderSettings] = None,
                display: Optional[DisplaySettings] = None, select=None,
                n_threads: int = 0, use_extension: bool = True) -> RenderedImage:
    """Render a localization table.

    ``select`` narrows the table to a subset: a :class:`~smappy.filter.LocFilter`,
    a boolean mask, or an index array.  Only the three or four columns the render
    needs are indexed, so filtering never copies the table.
    """
    settings = settings or RenderSettings()
    display = display or DisplaySettings()
    select = getattr(select, "indices", select)  # a LocFilter, or a mask/indices

    x, y = positions(locs, select)

    if settings.mode == "hist":
        sigma = None
    elif settings.mode == "gauss":
        sigma = settings.sigma
    elif settings.mode == "precision":
        name = _pick(locs, PRECISION_FIELDS, "localization precision")
        precision = locs[name] if select is None else locs[name][select]
        sigma = settings.sigma_settings.apply(precision, fov.pixelsize)
    else:
        raise ValueError(f"unknown mode {settings.mode!r}; "
                         "use 'hist', 'gauss' or 'precision'")

    weights = None
    if settings.weight_field is not None:
        weights = locs[settings.weight_field]
        weights = weights if select is None else weights[select]

    colors = None
    if settings.color_field is not None:
        values = locs[settings.color_field]
        values = values if select is None else values[select]
        lo, hi = settings.color_range or (float(np.nanmin(values)),
                                          float(np.nanmax(values)))
        colors = luts.colors(values, display.lut, lo, hi, display.invert)

    return render(x, y, fov, sigma=sigma, weights=weights, colors=colors,
                  roi_sigma=settings.roi_sigma, n_threads=n_threads,
                  use_extension=use_extension)
