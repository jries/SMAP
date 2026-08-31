"""The grid index must be complete: never miss a localization in the rectangle."""
import numpy as np
import pytest

from smappy.spatial import SpatialIndex


def _points(n=20000, seed=0):
    rng = np.random.default_rng(seed)
    # clustered, so cells are unevenly filled
    centers = rng.uniform(0, 1000, (20, 2))
    which = rng.integers(0, len(centers), n)
    return (centers[which] + rng.normal(0, 30, (n, 2))).T.astype(np.float32)


@pytest.mark.parametrize("seed", range(4))
def test_queries_are_complete(seed):
    x, y = _points()
    index = SpatialIndex(x, y)
    rng = np.random.default_rng(seed)
    for _ in range(20):
        x0, x1 = np.sort(rng.uniform(-100, 1100, 2))
        y0, y1 = np.sort(rng.uniform(-100, 1100, 2))
        found = index.query(x0, x1, y0, y1)
        truth = np.flatnonzero((x >= x0) & (x <= x1) & (y >= y0) & (y <= y1))
        assert np.isin(truth, found).all()          # nothing inside is missed
        assert len(np.unique(found)) == len(found)  # and nothing is duplicated


def test_the_margin_reaches_beyond_the_rectangle():
    x = np.array([0.0, 50.0, 100.0], np.float32)
    y = np.zeros(3, np.float32)
    index = SpatialIndex(x, y, cell_size=10.0)
    assert index.query(45, 55, -5, 5).size >= 1
    assert set(index.query(45, 55, -5, 5, margin=60).tolist()) == {0, 1, 2}


def test_non_finite_positions_are_never_returned():
    x, y = _points(n=1000, seed=1)
    x[::10] = np.nan
    y[5::10] = np.inf
    index = SpatialIndex(x, y)
    found = index.query(-1e6, 1e6, -1e6, 1e6)
    assert np.isfinite(x[found]).all() and np.isfinite(y[found]).all()
    assert found.size == int((np.isfinite(x) & np.isfinite(y)).sum())


def test_degenerate_inputs():
    empty = SpatialIndex(np.array([], np.float32), np.array([], np.float32))
    assert empty.query(-1, 1, -1, 1).size == 0
    same = SpatialIndex(np.full(100, 5.0, np.float32), np.full(100, 5.0, np.float32))
    assert same.query(4, 6, 4, 6).size == 100
    assert same.query(10, 20, 10, 20).size == 0
