"""Grouping: merging localizations of one emitter across consecutive frames.

A single emitter is localized in many consecutive frames, so a raw table counts
one blink many times.  Bright, long-lived emitters become hot spots in a
histogram render, and every statistic downstream sees a cluster of "independent"
points where there is one molecule.  Grouping collapses each such run into one
localization -- with a better precision, since it averages several measurements.

Two parts, following SMAP's ``Grouper.m`` and ``connectsingle2c.c``:

``connect``   greedy frame-to-frame linking (in C++, see ``csrc/group.hpp``).
``combine``   reduce each group to one row.

The combine rules are SMAP's, which are per column and not arbitrary:

===============================  ==========================================
positions, PSF size, colour      weighted mean, ``w = 1 / precision^2``
photons, background              sum
precisions and ``*_err``         ``1 / sqrt(sum(1 / e^2))``
log-likelihood                   max (the best frame's)
frame                            min (where the group starts)
===============================  ==========================================

Two departures, both because the blanket rule is wrong for the column:

* the error of a **summed** quantity adds in quadrature, ``sqrt(sum(e^2))``, not
  by the precision rule -- so ``photons_err`` and ``background_err`` use that.
  This is the rule that keeps shot noise consistent: with ``e_i = sqrt(N_i)`` it
  gives ``sqrt(sum(N_i))`` exactly, which is the shot noise of the summed
  photons.  SMAP applies its ``*_err`` rule to them instead, understating the
  error of a four-member group by a factor of six.
* each coordinate is weighted by **its own** error -- ``x_nm`` by ``x_err_nm``,
  ``y_nm`` by ``y_err_nm``, ``z_nm`` by ``z_err_nm`` -- where the table has
  them.  SMAP weights all three with the pooled lateral precision because that
  is the one weight it computes, and flags it as a shortcut; under astigmatism
  ``x_err`` and ``y_err`` diverge with z (that divergence is what encodes z), so
  the pooled weight is the right one only when they happen to be equal.
* ``logl_rel`` is a likelihood *per pixel* and stays comparable under max, but
  the raw ``logl`` of a group is the sum over its members' fits; taking the max
  of it is meaningless across groups of different size, so it is dropped rather
  than given a wrong value.  ``logl_rel`` -- which is what the filter uses --
  is kept.

Grouping is expensive (the linking is sequential and cannot be vectorised), so
it is done once and the grouped table kept alongside the original rather than
recomputed when the display switches between them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .locs import Localizations

try:
    from . import _group
except ImportError:  # pure-Python fallback if the extension is not built
    _group = None


# How each column is combined.  Taken from SMAP's `Grouper.m`; names translated
# to this package's, and the two corrections described in the module docstring.
COMBINE_MODES: Dict[str, str] = {
    "x_nm": "mean", "y_nm": "mean", "z_nm": "mean",
    "x_pix": "mean", "y_pix": "mean",
    "peak_x_nm": "mean", "peak_y_nm": "mean",
    "peak_x_pix": "mean", "peak_y_pix": "mean",
    "sigma_nm": "mean", "sigma_pix": "mean",
    "sigma_x_nm": "mean", "sigma_y_nm": "mean",
    "sigma_x_pix": "mean", "sigma_y_pix": "mean",
    "photons": "sum", "background": "sum",
    "photons_err": "quad", "background_err": "quad",
    "loc_precision_nm": "precision", "loc_precision_pix": "precision",
    "x_err_nm": "precision", "y_err_nm": "precision", "z_err_nm": "precision",
    "x_err_pix": "precision", "y_err_pix": "precision",
    "logl_rel": "max",
    "frame": "min",
    "iterations": "max",
}

# Columns that have no meaningful group value; see the module docstring.
DROP_ON_GROUPING = ("logl",)

# The precision columns the general weight is taken from, best first (SMAP's
# order).  It is a lateral precision, so it is the right weight for x and y.
WEIGHT_FIELDS = ("loc_precision_nm", "loc_precision_pix", "x_err_nm", "x_err_pix")

# Columns that carry their own uncertainty and are weighted by it rather than by
# the pooled lateral precision.  SMAP weights everything with `locprecnm`
# because that is the one weight it computes -- a comment in `Grouper.m` flags
# it as a shortcut -- but the fitter returns a separate error per coordinate and
# under astigmatism they genuinely differ: `x_err` and `y_err` diverge with z
# (that is what encodes z), and `z_err` varies by a factor of several over the
# range.  Weighting each coordinate by its own error is the correct inverse-
# variance estimate; the pooled one is only right when the errors are equal.
WEIGHT_FOR: Dict[str, Tuple[str, ...]] = {
    "x_nm": ("x_err_nm",), "y_nm": ("y_err_nm",), "z_nm": ("z_err_nm",),
    "x_pix": ("x_err_pix",), "y_pix": ("y_err_pix",),
}


def _mode(name: str) -> str:
    if name in COMBINE_MODES:
        return COMBINE_MODES[name]
    return "precision" if name.endswith(("_err", "err")) else "mean"


def connect(x, y, frame, dx: float = 50.0, dt: int = 1,
            blocks: Optional[np.ndarray] = None) -> np.ndarray:
    """Assign every localization a 1-based group id, in the input order.

    ``dx`` is the half-width of the search box in the units of x and y, ``dt``
    the number of dark frames a particle may skip.  ``blocks`` labels groups of
    localizations that linking may not cross (a file or channel number); linking
    is run once per block, rather than SMAP's trick of zeroing the frame at each
    boundary, which leaves the array no longer sorted by frame.
    """
    if _group is None:
        raise RuntimeError("the _group extension is not built; "
                           "run setup.py build_ext --build-lib src")
    x = np.asarray(x, np.float64)
    y = np.asarray(y, np.float64)
    frame = np.asarray(np.rint(np.asarray(frame, np.float64)), np.int64)
    if not (x.shape == y.shape == frame.shape) or x.ndim != 1:
        raise ValueError("x, y and frame must be 1-D and the same length")

    keys: Tuple[np.ndarray, ...] = ()
    if blocks is not None:
        b = np.asarray(blocks)
        keys = tuple(b.T) if b.ndim == 2 else (b,)

    # the linker needs (frame, x) ascending, within a block
    order = np.lexsort((x, frame) + tuple(reversed(keys)))
    out = np.zeros(x.size, np.int64)

    if keys:
        stacked = np.stack([np.asarray(k)[order] for k in keys], axis=1)
        edges = np.concatenate(
            ([0], np.flatnonzero(np.any(stacked[1:] != stacked[:-1], axis=1)) + 1,
             [x.size]))
    else:
        edges = np.array([0, x.size])

    offset = 0
    for begin, end in zip(edges[:-1], edges[1:]):
        block = order[begin:end]
        ids, n_groups = _group.connect(x[block], y[block], frame[block], dx, dt)
        out[block] = ids + offset
        offset += n_groups
    return out


def _inverse_variance(values) -> np.ndarray:
    w = 1.0 / np.asarray(values, np.float64) ** 2
    return np.where(np.isfinite(w), w, 1.0)           # SMAP: infinite weight -> 1


def _weights(locs: Localizations) -> np.ndarray:
    """Per-localization weight for the means: ``1 / precision^2``, SMAP's order."""
    for name in WEIGHT_FIELDS:
        if name in locs:
            return _inverse_variance(locs[name])
    if "photons" in locs:
        return np.asarray(locs["photons"], np.float64)
    return np.ones(len(locs))


def combine(locs: Localizations, group_index: np.ndarray,
            fields: Optional[Sequence[str]] = None) -> Localizations:
    """Reduce each group to one row, one column at a time by its own rule."""
    gi = np.asarray(group_index, np.int64)
    if gi.size and gi.min() < 1:
        raise ValueError("group ids must be 1-based")
    n_groups = int(gi.max()) if gi.size else 0
    size = n_groups + 1

    default_w = _weights(locs)
    weights: Dict[str, np.ndarray] = {"": default_w}
    sums: Dict[str, np.ndarray] = {
        "": np.bincount(gi, weights=default_w, minlength=size)[1:]}
    n_in_group = np.bincount(gi, minlength=size)[1:]

    def weighting(name: str) -> Tuple[np.ndarray, np.ndarray]:
        """The weight for this column: its own uncertainty if it has one."""
        for column in WEIGHT_FOR.get(name, ()):
            if column in locs:
                if column not in weights:
                    weights[column] = _inverse_variance(locs[column])
                    sums[column] = np.bincount(gi, weights=weights[column],
                                               minlength=size)[1:]
                return weights[column], sums[column]
        return weights[""], sums[""]

    names = list(fields) if fields is not None else \
        [n for n in locs.keys() if n not in DROP_ON_GROUPING]

    # min and max need a per-group reduction; sorting once beats np.minimum.at,
    # which is unbuffered and slow
    extremes = any(_mode(n) in ("min", "max") for n in names)
    if extremes:
        order = np.argsort(gi, kind="stable")
        starts = np.concatenate(
            ([0], np.flatnonzero(gi[order][1:] != gi[order][:-1]) + 1))

    columns: Dict[str, np.ndarray] = {}
    for name in names:
        values = np.asarray(locs[name], np.float64)
        mode = _mode(name)
        if mode == "mean":
            w, sum_w = weighting(name)
            result = np.bincount(gi, weights=values * w, minlength=size)[1:] / sum_w
        elif mode == "sum":
            result = np.bincount(gi, weights=values, minlength=size)[1:]
        elif mode == "precision":
            result = 1.0 / np.sqrt(
                np.bincount(gi, weights=1.0 / values ** 2, minlength=size)[1:])
        elif mode == "quad":   # the error of a sum, not of a mean
            result = np.sqrt(
                np.bincount(gi, weights=values ** 2, minlength=size)[1:])
        elif mode in ("min", "max"):
            reduce = np.minimum if mode == "min" else np.maximum
            result = reduce.reduceat(values[order], starts)
        else:
            raise ValueError(f"unknown combine mode {mode!r} for {name!r}")
        columns[name] = result.astype(np.float32)

    for name in ("frame", "iterations"):
        if name in columns:
            columns[name] = columns[name].astype(np.int64 if name == "frame"
                                                 else np.int32)
    columns["n_in_group"] = n_in_group.astype(np.int32)

    metadata = dict(locs.metadata)
    metadata["grouped"] = True
    return Localizations(columns, metadata)


@dataclass(frozen=True)
class GroupSettings:
    """How to link.  ``dx`` is in the units of the table (nm, normally)."""

    dx: float = 50.0
    dt: int = 1
    block_fields: Sequence[str] = ("filenumber", "channel")


def group(locs: Localizations, settings: Optional[GroupSettings] = None
          ) -> Tuple[Localizations, np.ndarray]:
    """Group a table.  Returns the grouped table and the per-input group id.

    The ids are **1-based**, as `connect` produces them, so the row of the
    grouped table a localization ended up in is ``group_index - 1``.
    """
    settings = settings or GroupSettings()
    from .render import positions

    x, y = positions(locs)
    if "frame" not in locs:
        raise KeyError("grouping needs a 'frame' column")

    present = [locs[name] for name in settings.block_fields if name in locs]
    blocks = np.stack(present, axis=1) if present else None
    group_index = connect(x, y, locs["frame"], settings.dx, settings.dt, blocks)
    return combine(locs, group_index), group_index
