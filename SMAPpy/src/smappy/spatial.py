"""A grid index over localization positions.

Zooming in is the common case in a viewer, and without an index every render
walks all five million localizations to find the fifty thousand on screen --
which, once the rendering itself is threaded, is the dominant cost.  The index
answers "which localizations are in this rectangle" in time proportional to the
answer.

The structure is the simplest one that does that: localizations are sorted by
grid cell, and each cell's entries are a contiguous run.  Cell ids run along x
first, so a whole row of cells is *one* slice of the sorted order and a
rectangular query is one slice per cell row.  Building it is a single argsort;
it depends only on the positions, so it survives every change of filter,
settings and field of view.

The index is not a filter.  It returns candidates -- possibly a few outside the
rectangle, since a query is padded so blobs whose centre lies just outside still
contribute their tails -- and the renderer clips.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

# Cells per axis.  The cap keeps the sort keys 16-bit (see above); it is far
# above what the target occupancy asks for at any realistic size.
MAX_CELLS = 4096


class SpatialIndex:
    """Localizations bucketed into a uniform grid, sorted by cell."""

    def __init__(self, x, y, target_per_cell: int = 256,
                 cell_size: Optional[float] = None,
                 extent: Optional[Sequence[float]] = None):
        """``extent`` is ``(x0, y0, x1, y1)`` the grid must cover *at least*.

        The grid otherwise follows the data, which is right for a table that is
        complete but moves under a table that is still growing: pass the camera
        field of view and the index keeps the same bounds however little of it
        has been filled so far.  A larger extent is never wrong, only coarser.
        """
        self.x = np.ascontiguousarray(x, dtype=np.float32)
        self.y = np.ascontiguousarray(y, dtype=np.float32)
        n = self.x.size

        finite = np.isfinite(self.x) & np.isfinite(self.y)
        if not finite.any() and extent is None:
            self.x0 = self.y0 = 0.0
            self.cell_size = 1.0
            self.n_cols = self.n_rows = 1
            self.order = np.empty(0, np.int64)
            self.start = np.zeros(2, np.int64)
            self._end = 0
            self._bounds = (0.0, 0.0, 1.0, 1.0)
            return

        if finite.any():
            x0, x1 = float(self.x[finite].min()), float(self.x[finite].max())
            y0, y1 = float(self.y[finite].min()), float(self.y[finite].max())
        else:
            x0 = y0 = np.inf
            x1 = y1 = -np.inf
        if extent is not None:
            ex0, ey0, ex1, ey1 = (float(v) for v in extent)
            x0, y0 = min(x0, ex0), min(y0, ey0)
            x1, y1 = max(x1, ex1), max(y1, ey1)

        self.x0, self.y0 = x0, y0
        # what the index actually covers, which is not the grid: the grid is
        # rounded up to a whole cell, and one cell can be enormous when a
        # block holds a handful of localizations
        self._bounds = (x0, y0, x1, y1)
        width, height = x1 - x0, y1 - y0

        if cell_size is None:
            # aim for `target_per_cell` localizations in a cell, square in data
            # space; the point is to keep both the cell count and the number of
            # candidates per query moderate
            area = max(width * height, np.finfo(np.float32).tiny)
            cell_size = float(np.sqrt(area * target_per_cell / max(n, 1)))
        self.cell_size = max(cell_size, np.finfo(np.float32).tiny,
                             max(width, height) / MAX_CELLS)

        self.n_cols = max(int(width / self.cell_size) + 1, 1)
        self.n_rows = max(int(height / self.cell_size) + 1, 1)

        col = self._cell_key(self.x, self.x0, self.n_cols, finite)
        row = self._cell_key(self.y, self.y0, self.n_rows, finite)
        # non-finite positions get a sentinel row past the end of the grid, so
        # they sort last and no query can reach them: they cannot be rendered
        np.copyto(row, self.n_rows, where=~finite)

        self.order = np.lexsort((col, row))
        n_cells = self.n_cols * self.n_rows
        counts = np.bincount(row.astype(np.int64) * self.n_cols + col,
                             minlength=n_cells + 1)
        self.start = np.concatenate([[0], np.cumsum(counts[:n_cells])])
        self._end = int(self.start[n_cells])   # where the non-finite tail begins

    def _cell_key(self, values, origin, n_cells, finite) -> np.ndarray:
        cells = (values - origin) * (1.0 / self.cell_size)
        return np.clip(np.where(finite, cells, 0), 0, n_cells - 1).astype(np.uint16)

    @property
    def n_localizations(self) -> int:
        return self.x.size

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The rectangle covered: the localizations, and any ``extent`` given."""
        return self._bounds

    def query(self, x0: float, x1: float, y0: float, y1: float,
              margin: float = 0.0) -> np.ndarray:
        """Indices of the localizations in the rectangle, plus a few outside."""
        c0 = int(np.floor((x0 - margin - self.x0) / self.cell_size))
        c1 = int(np.floor((x1 + margin - self.x0) / self.cell_size))
        r0 = int(np.floor((y0 - margin - self.y0) / self.cell_size))
        r1 = int(np.floor((y1 + margin - self.y0) / self.cell_size))
        c0, c1 = max(c0, 0), min(c1, self.n_cols - 1)
        r0, r1 = max(r0, 0), min(r1, self.n_rows - 1)
        if c0 > c1 or r0 > r1:
            return np.empty(0, np.int64)

        # the whole grid: hand back everything rather than stitching every row
        if c0 == 0 and c1 == self.n_cols - 1 and r0 == 0 and r1 == self.n_rows - 1:
            return self.order[:self._end]

        rows = np.arange(r0, r1 + 1) * self.n_cols
        begins = self.start[rows + c0]
        ends = self.start[rows + c1 + 1]
        if len(rows) == 1:
            return self.order[begins[0]:ends[0]]
        return np.concatenate([self.order[b:e] for b, e in zip(begins, ends)
                               if e > b]) if (ends > begins).any() \
            else np.empty(0, np.int64)


class GrowingIndex:
    """A spatial index over a table that is still being appended to.

    Rebuilding from scratch on every update would dominate it: the argsort over
    five million localizations costs a third of a second, and an online fit
    delivers a block every few seconds.  So each appended block gets its own
    index, a query is the union of what the segments answer, and segments are
    merged only occasionally -- the tail into one when there are too many of
    them, everything into one once the tail has grown to the size of the head.
    That keeps the total merge work O(N log N) over an acquisition instead of
    O(N) per block, while a query stays a handful of slices per segment.

    Segments deliberately do *not* share a grid.  Each `query` is answered
    against the segment's own origin and cell size, so a block that lands
    outside the area seen so far, or one far sparser than the last, is indexed
    on its own terms and never forces a rebuild.  Passing ``extent`` (the
    camera's field of view) still helps: it pins `bounds`, and with it the
    viewer's full view, so the image does not rescale as data arrives.

    Indices are positions in the appended-to table, so they stay valid as it
    grows -- which is why blocks may only be added at the end.
    """

    def __init__(self, x, y, extent: Optional[Sequence[float]] = None,
                 target_per_cell: int = 256, max_segments: int = 8):
        self.extent = None if extent is None else tuple(float(v) for v in extent)
        self.target_per_cell = int(target_per_cell)
        self.max_segments = max(int(max_segments), 1)
        self._x: List[np.ndarray] = []
        self._y: List[np.ndarray] = []
        self._segments: List[SpatialIndex] = []
        self._offsets: List[int] = []
        self._n = 0
        self.append(x, y)

    # ------------------------------------------------------------------ growing
    def append(self, x, y) -> "GrowingIndex":
        """Index another block, whose rows follow the ones already indexed."""
        x = np.ascontiguousarray(x, dtype=np.float32)
        y = np.ascontiguousarray(y, dtype=np.float32)
        if x.size != y.size:
            raise ValueError(f"{x.size} x against {y.size} y")
        if x.size == 0:
            return self

        self._offsets.append(self._n)
        self._x.append(x)
        self._y.append(y)
        self._segments.append(SpatialIndex(x, y, self.target_per_cell,
                                           extent=self.extent))
        self._n += x.size

        if len(self._segments) > self.max_segments:
            head = self._segments[0].n_localizations
            self._merge(0 if self._n - head >= head else 1)
        return self

    def _merge(self, start: int) -> None:
        """Replace segments ``start:`` by one index over the same rows."""
        if len(self._segments) - start < 2:
            return
        x = np.concatenate(self._x[start:])
        y = np.concatenate(self._y[start:])
        self._x[start:] = [x]
        self._y[start:] = [y]
        self._offsets[start:] = [self._offsets[start]]
        self._segments[start:] = [SpatialIndex(x, y, self.target_per_cell,
                                               extent=self.extent)]

    # ------------------------------------------------------------------ reading
    @property
    def n_localizations(self) -> int:
        return self._n

    @property
    def n_segments(self) -> int:
        return len(self._segments)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The rectangle every segment covers together."""
        boxes = [s.bounds for s in self._segments]
        if self.extent is not None:
            boxes.append(self.extent)
        if not boxes:
            return (0.0, 0.0, 1.0, 1.0)
        corners = np.array(boxes, dtype=np.float64)
        return (float(corners[:, 0].min()), float(corners[:, 1].min()),
                float(corners[:, 2].max()), float(corners[:, 3].max()))

    def query(self, x0: float, x1: float, y0: float, y1: float,
              margin: float = 0.0) -> np.ndarray:
        """Indices of the localizations in the rectangle, plus a few outside."""
        parts = []
        for segment, offset in zip(self._segments, self._offsets):
            found = segment.query(x0, x1, y0, y1, margin)
            if found.size:
                parts.append(found + offset if offset else found)
        if not parts:
            return np.empty(0, np.int64)
        return parts[0] if len(parts) == 1 else np.concatenate(parts)
