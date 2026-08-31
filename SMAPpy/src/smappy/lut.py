"""Look-up tables: 256-entry colour ramps.

A LUT is an ``(256, 3)`` float32 array of RGB in [0, 1].  Two places use them
and they are not the same operation:

* an intensity image is *recoloured* -- one LUT lookup per output pixel, after
  normalisation, so changing the LUT costs nothing and never touches the
  localizations;
* a field-coloured image needs the colour of every localization *before*
  accumulation, because two localizations in one pixel must mix.

Most ramps are defined by a handful of anchor colours interpolated to 256
entries, which is what SMAP's ``usercolormap`` does.  ``hot``, ``jet``, ``hsv``
and ``gray`` follow MATLAB's definitions so images look the same as SMAP's.
No matplotlib dependency: the tables are small and the interpolation is three
lines.
"""

from __future__ import annotations

from typing import Dict, Sequence, Union

import numpy as np

N = 256  # entries in every LUT


def from_anchors(*colors: Sequence[float]) -> np.ndarray:
    """A ramp through equally spaced anchor colours (SMAP's ``usercolormap``)."""
    anchors = np.asarray(colors, dtype=np.float64)
    at = np.linspace(0.0, 1.0, len(anchors))
    x = np.linspace(0.0, 1.0, N)
    return np.stack([np.interp(x, at, anchors[:, k]) for k in range(3)],
                    axis=1).astype(np.float32)


def _matlab_hot() -> np.ndarray:
    n = int(np.ceil(N * 3 / 8))
    ramp = np.arange(1, n + 1) / n
    r = np.concatenate([ramp, np.ones(N - n)])
    g = np.concatenate([np.zeros(n), ramp, np.ones(N - 2 * n)])
    b = np.concatenate([np.zeros(2 * n), np.arange(1, N - 2 * n + 1) / (N - 2 * n)])
    return np.stack([r, g, b], axis=1).astype(np.float32)


def _matlab_jet() -> np.ndarray:
    x = np.linspace(0.0, 1.0, N)
    def band(center):  # trapezoid of half-width 1/8, plateau 1/8 each side
        return np.clip(1.5 - np.abs(x - center) * 4.0, 0.0, 1.0)
    return np.stack([band(0.6875), band(0.5), band(0.3125)],
                    axis=1).astype(np.float32)


def _matlab_hsv() -> np.ndarray:
    h = np.linspace(0.0, 1.0, N, endpoint=False) * 6.0
    i, f = np.floor(h).astype(int) % 6, h - np.floor(h)
    up, down = f, 1.0 - f
    one, zero = np.ones_like(f), np.zeros_like(f)
    r = np.choose(i, [one, down, zero, zero, up, one])
    g = np.choose(i, [up, one, one, down, zero, zero])
    b = np.choose(i, [zero, zero, up, one, one, down])
    return np.stack([r, g, b], axis=1).astype(np.float32)


def _single(r: float, g: float, b: float) -> np.ndarray:
    ramp = np.linspace(0.0, 1.0, N)
    return np.stack([ramp * r, ramp * g, ramp * b], axis=1).astype(np.float32)


_gray = _single(1, 1, 1)

# Turbo, sampled from SMAP's copy of the table every 8th entry.  It is defined
# by a table rather than a formula, and 33 anchors reproduce it to 3/255 -- less
# than a step of the 256-entry ramp -- so it is stored as anchors like the rest.
# Worth having for z: unlike `hsv` it does not wrap, so the two ends of the
# range cannot be confused, and unlike `jet` its lightness increases
# monotonically, so structure does not appear where the colour map bunches up.
_TURBO_ANCHORS = [
    [0.1900, 0.0718, 0.2322],
    [0.2250, 0.1635, 0.4510],
    [0.2511, 0.2524, 0.6337],
    [0.2682, 0.3382, 0.7805],
    [0.2763, 0.4212, 0.8912],
    [0.2754, 0.5011, 0.9659],
    [0.2586, 0.5796, 0.9988],
    [0.2138, 0.6589, 0.9796],
    [0.1584, 0.7355, 0.9231],
    [0.1117, 0.8057, 0.8452],
    [0.0927, 0.8655, 0.7623],
    [0.1201, 0.9119, 0.6866],
    [0.1966, 0.9490, 0.5947],
    [0.3051, 0.9770, 0.4899],
    [0.4278, 0.9942, 0.3857],
    [0.5466, 0.9991, 0.2958],
    [0.6436, 0.9900, 0.2336],
    [0.7260, 0.9647, 0.2064],
    [0.8047, 0.9245, 0.2046],
    [0.8753, 0.8727, 0.2155],
    [0.9330, 0.8124, 0.2267],
    [0.9732, 0.7468, 0.2254],
    [0.9931, 0.6741, 0.2035],
    [0.9959, 0.5870, 0.1690],
    [0.9836, 0.4929, 0.1285],
    [0.9580, 0.3996, 0.0883],
    [0.9211, 0.3149, 0.0548],
    [0.8742, 0.2453, 0.0330],
    [0.8161, 0.1846, 0.0181],
    [0.7462, 0.1310, 0.0085],
    [0.6645, 0.0844, 0.0042],
    [0.5710, 0.0447, 0.0053],
    [0.4796, 0.0158, 0.0106],
]

# The subset of SMAP's LUTs that is actually used: a few intensity ramps for
# recolouring, and a few full-hue ramps for field colouring.
_TABLES: Dict[str, np.ndarray] = {
    "hot": _matlab_hot(),
    "gray": _gray,
    "gray_inverted": _gray[::-1].copy(),
    "jet": _matlab_jet(),
    "turbo": from_anchors(*_TURBO_ANCHORS),
    "hsv": _matlab_hsv(),
    "red": _single(1, 0, 0),
    "green": _single(0, 1, 0),
    "blue": _single(0, 0, 1),
    "cyan": _single(0, 1, 1),
    "magenta": _single(1, 0, 1),
    "yellow": _single(1, 1, 0),
    "green_cold": from_anchors([0, 0, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1]),
    "cyan_cold": from_anchors([0, 0, 0], [0, .2, 1], [0, .5, 1], [0, 1, 1],
                              [.5, 1, 1], [1, 1, 1]),
    "kgy": from_anchors([0, 0, 0], [.5, .25, 0], [1, .5, 0], [1, 1, 0], [1, 1, 1]),
    "bgy": from_anchors([0, 0, 1], [0, 1, 0], [1, 1, 0]),
    "bry": from_anchors([0, 0, 1], [1, 0, 0], [1, 1, 0]),
}

LUT = Union[str, np.ndarray]


def names() -> list:
    return sorted(_TABLES)


def get(lut: LUT, invert: bool = False) -> np.ndarray:
    """Resolve a LUT name (or pass an array through), optionally reversed."""
    if isinstance(lut, str):
        try:
            table = _TABLES[lut]
        except KeyError:
            raise KeyError(f"unknown LUT {lut!r}; have {', '.join(names())}") from None
    else:
        table = np.asarray(lut, dtype=np.float32)
        if table.ndim != 2 or table.shape[1] != 3:
            raise ValueError("a LUT must have shape (n, 3)")
    return table[::-1].copy() if invert else table


def register(name: str, table: np.ndarray) -> None:
    """Add a LUT, e.g. one taken from matplotlib."""
    _TABLES[name] = get(table)


def colors(values: np.ndarray, lut: LUT, vmin: float, vmax: float,
           invert: bool = False) -> np.ndarray:
    """Map field values to ``(n, 3)`` colours, clamped to the range."""
    table = get(lut, invert)
    scale = (len(table) - 1) / (vmax - vmin) if vmax > vmin else 0.0
    index = np.clip((np.asarray(values, dtype=np.float32) - vmin) * scale,
                    0, len(table) - 1).astype(np.int32)
    return table[index]
