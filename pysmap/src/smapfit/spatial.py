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

from typing import Optional

import numpy as np

# Cells per axis.  The cap keeps the sort keys 16-bit (see above); it is far
# above what the target occupancy asks for at any realistic size.
MAX_CELLS = 4096


class SpatialIndex:
    """Localizations bucketed into a uniform grid, sorted by cell."""

    def __init__(self, x, y, target_per_cell: int = 256,
                 cell_size: Optional[float] = None):
        self.x = np.ascontiguousarray(x, dtype=np.float32)
        self.y = np.ascontiguousarray(y, dtype=np.float32)
        n = self.x.size

        finite = np.isfinite(self.x) & np.isfinite(self.y)
        if not finite.any():
            self.x0 = self.y0 = 0.0
            self.cell_size = 1.0
            self.n_cols = self.n_rows = 1
            self.order = np.empty(0, np.int64)
            self.start = np.zeros(2, np.int64)
            self._end = 0
            return

        self.x0 = float(self.x[finite].min())
        self.y0 = float(self.y[finite].min())
        width = float(self.x[finite].max()) - self.x0
        height = float(self.y[finite].max()) - self.y0

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
