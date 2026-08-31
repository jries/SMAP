"""Selecting which localizations to display.

Filtering is a range per field -- localization precision below 30 nm, relative
log-likelihood above -2, z within a slab -- and the displayed set is the
intersection of all of them.  Recomputing that from scratch on every render is
what makes interactive filtering slow, so, as in SMAP, **one boolean array is
cached per field** and only the field that changed is recomputed; the
intersection is then an AND over a handful of arrays, which is memory-bandwidth
work and takes a few milliseconds even for millions of localizations.

The result is exposed two ways.  ``mask`` is the boolean array, for combining
with other selections; ``indices`` is the compacted version, which is what the
renderer wants -- it indexes three or four columns instead of copying the whole
table.

NaN is excluded by any range, because a comparison against NaN is false.  That
is deliberate: a localization with no z should not appear in a z slab.

Not ported from SMAP: ROI/polygon masks and the grouped-localization filters.
``set_mask`` takes an arbitrary boolean array under a name, which is where such
a selection -- or the viewer's "inside the field of view" -- would plug in.
"""

from __future__ import annotations

from typing import Dict, Iterator, Optional, Tuple

import numpy as np

from .locs import Localizations


class LocFilter:
    """Cached range filters over the fields of a :class:`Localizations` table."""

    def __init__(self, locs: Localizations, **ranges: Tuple[float, float]):
        self._locs = locs
        self._n = len(locs)
        self._ranges: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
        self._masks: Dict[str, np.ndarray] = {}
        self._combined: Optional[np.ndarray] = None
        self._indices: Optional[np.ndarray] = None
        for field, (lo, hi) in ranges.items():
            self.set(field, lo, hi)

    # ------------------------------------------------------------- setting up
    def set(self, field: str, lo: Optional[float] = None,
            hi: Optional[float] = None) -> "LocFilter":
        """Keep localizations with ``lo <= field <= hi``; ``None`` is unbounded.

        Only this field's array is recomputed.
        """
        self._ranges[field] = (lo, hi)
        self._masks[field] = self._range_mask(field, lo, hi)
        return self._invalidate()

    def _range_mask(self, field: str, lo: Optional[float], hi: Optional[float],
                    start: int = 0) -> np.ndarray:
        """The mask for one range, over the rows from ``start`` on."""
        values = self._locs[field][start:]  # raises if the column does not exist
        mask = np.ones(values.shape[0], dtype=bool)
        if lo is not None:
            mask &= values >= lo
        if hi is not None:
            mask &= values <= hi
        return mask

    def append(self, n_new: int) -> "LocFilter":
        """Take in ``n_new`` rows appended to the table since the last call.

        Only the new rows are tested, against the bounds already set: appending
        never re-derives a limit, so what the user typed keeps meaning exactly
        what they typed while data arrives.  A mask set by hand (`set_mask`)
        has no rule to extend it, so its new rows are excluded until it is set
        again -- visible and recoverable, unlike silently including them.
        """
        start, self._n = self._n, self._n + int(n_new)
        if len(self._locs) != self._n:
            raise ValueError(f"the table has {len(self._locs)} rows, "
                             f"expected {self._n} after appending {n_new}")
        for name, mask in self._masks.items():
            tail = (self._range_mask(name, *self._ranges[name], start=start)
                    if name in self._ranges
                    else np.zeros(int(n_new), dtype=bool))
            self._masks[name] = np.concatenate([mask, tail])
        return self._invalidate()

    def set_mask(self, name: str, mask: np.ndarray) -> "LocFilter":
        """Add an arbitrary selection under a name (a ROI, a manual pick)."""
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (self._n,):
            raise ValueError(f"a mask must have {self._n} entries, got {mask.shape}")
        self._ranges.pop(name, None)
        self._masks[name] = mask
        return self._invalidate()

    def remove(self, name: str) -> "LocFilter":
        self._masks.pop(name, None)
        self._ranges.pop(name, None)
        return self._invalidate()

    def clear(self) -> "LocFilter":
        self._masks.clear()
        self._ranges.clear()
        return self._invalidate()

    def _invalidate(self) -> "LocFilter":
        self._combined = None
        self._indices = None
        return self

    # --------------------------------------------------------------- the result
    @property
    def mask(self) -> np.ndarray:
        """Boolean array over the table: the intersection of every filter."""
        if self._combined is None:
            if not self._masks:
                self._combined = np.ones(self._n, dtype=bool)
            else:
                self._combined = np.logical_and.reduce(list(self._masks.values()))
        return self._combined

    @property
    def indices(self) -> np.ndarray:
        """The positions of the kept localizations; what the renderer indexes with."""
        if self._indices is None:
            self._indices = np.flatnonzero(self.mask)
        return self._indices

    @property
    def ranges(self) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
        return dict(self._ranges)

    def __len__(self) -> int:
        return int(self.mask.sum())

    def __iter__(self) -> Iterator[str]:
        return iter(self._masks)

    def __contains__(self, name: str) -> bool:
        return name in self._masks

    def apply(self) -> Localizations:
        """The filtered table.  Copies every column; prefer ``indices``."""
        return self._locs[self.mask]

    def counts(self) -> Dict[str, int]:
        """How many localizations each filter keeps on its own.

        What a GUI shows next to each slider, and the quickest way to see which
        filter is doing the cutting.
        """
        return {name: int(mask.sum()) for name, mask in self._masks.items()}

    def __str__(self) -> str:
        if not self._masks:
            return f"no filter, {self._n} localizations"
        parts = []
        for name in self._masks:
            lo, hi = self._ranges.get(name, (None, None))
            if name in self._ranges:
                parts.append(f"{name} "
                             f"[{'' if lo is None else f'{lo:g}'}"
                             f"..{'' if hi is None else f'{hi:g}'}]")
            else:
                parts.append(f"{name} (mask)")
        return (f"{len(self)} of {self._n} localizations: " + ", ".join(parts))


def quantile_range(locs: Localizations, field: str, lower: float = 0.001,
                   upper: float = 0.999) -> Tuple[float, float]:
    """A sensible starting range for a field, ignoring its tails.

    Used to initialise filter limits and colour ranges without a stray outlier
    setting the scale.
    """
    values = np.asarray(locs[field], dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (0.0, 0.0)
    return (float(np.quantile(values, lower)), float(np.quantile(values, upper)))
