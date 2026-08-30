"""Appending to a live view must change nothing except the data.

Every property tested here is the same one: what you get after appending in
blocks is what you would have got from the whole table at once -- the same
selection, the same image -- while the settings the user made stay untouched.
"""
import numpy as np
import pytest

from smapfit.filter import LocFilter
from smapfit.locs import Localizations
from smapfit.render import FieldOfView, RenderSettings, render_locs
from smapfit.spatial import GrowingIndex, SpatialIndex
from smapfit.viewer import ViewState


def _table(n=20000, seed=0, first_frame=0):
    rng = np.random.default_rng(seed)
    return Localizations({
        "x_nm": rng.uniform(0, 10000, n).astype(np.float32),
        "y_nm": rng.uniform(0, 8000, n).astype(np.float32),
        "z_nm": rng.uniform(-400, 400, n).astype(np.float32),
        "loc_precision_nm": (rng.gamma(4, 3, n) + 5).astype(np.float32),
        "logl_rel": rng.normal(-1, 0.5, n).astype(np.float32),
        "frame": np.arange(first_frame, first_frame + n, dtype=np.int64),
    }, {"units": "nm"})


def _blocks(locs, sizes):
    start = 0
    for size in sizes:
        yield locs[start:start + size]
        start += size


# ------------------------------------------------------------------ the table
def test_extending_in_blocks_gives_the_same_table():
    whole = _table(5000)
    grown = Localizations()
    for block in _blocks(whole, [1, 700, 0, 2299, 2000]):
        grown.extend(block)
    assert len(grown) == len(whole)
    for name in whole.keys():
        assert grown[name].dtype == whole[name].dtype
        assert np.array_equal(grown[name], whole[name])


def test_a_column_is_never_longer_than_the_table():
    """The spare capacity must not be visible: columns are exact-length views."""
    grown = Localizations()
    for block in _blocks(_table(1000), [100] * 10):
        grown.extend(block)
        assert all(len(v) == len(grown) for v in grown.columns.values())
    grown.compact()
    assert all(v.base is None for v in grown.columns.values())


def test_extending_promotes_the_dtype_rather_than_truncating():
    grown = Localizations({"a": np.arange(3, dtype=np.int32)})
    grown.extend(Localizations({"a": np.array([0.5, 1.5])}))
    assert np.allclose(grown["a"], [0, 1, 2, 0.5, 1.5])


# ------------------------------------------------------------------ the index
@pytest.mark.parametrize("seed", range(3))
def test_the_growing_index_answers_what_a_rebuilt_one_would(seed):
    locs = _table(30000, seed=seed)
    x, y = np.asarray(locs["x_nm"]), np.asarray(locs["y_nm"])
    index = GrowingIndex(x[:100], y[:100], max_segments=3)
    for start in range(100, len(x), 1500):
        index.append(x[start:start + 1500], y[start:start + 1500])
    assert index.n_localizations == len(x)
    assert index.n_segments <= 4          # merged, not one segment per block

    rng = np.random.default_rng(seed)
    for _ in range(20):
        x0, x1 = np.sort(rng.uniform(-1000, 11000, 2))
        y0, y1 = np.sort(rng.uniform(-1000, 9000, 2))
        found = index.query(x0, x1, y0, y1)
        truth = np.flatnonzero((x >= x0) & (x <= x1) & (y >= y0) & (y <= y1))
        assert np.isin(truth, found).all()          # nothing inside is missed
        assert len(np.unique(found)) == len(found)  # and nothing is duplicated


def test_a_fixed_extent_keeps_the_bounds_still():
    """The whole point of passing the camera field of view: the view holds."""
    extent = (0.0, 0.0, 10000.0, 8000.0)
    locs = _table(2000)
    index = GrowingIndex(locs["x_nm"][:10], locs["y_nm"][:10], extent=extent)
    before = index.bounds
    index.append(locs["x_nm"][10:], locs["y_nm"][10:])
    assert index.bounds == before
    assert before[0] <= 0.0 and before[2] >= 10000.0


def test_without_an_extent_the_index_still_covers_what_arrives_late():
    index = GrowingIndex([0.0], [0.0])
    index.append([1000.0], [1000.0])
    assert index.query(990, 1010, 990, 1010).tolist() == [1]


# ----------------------------------------------------------------- the filter
def test_appending_extends_the_filter_without_rederiving_it():
    whole = _table(6000)
    part = whole[:2000]
    grown = Localizations()
    grown.extend(part)
    filt = LocFilter(grown, loc_precision_nm=(None, 15.0), z_nm=(-100.0, 100.0))
    ranges = filt.ranges

    for block in _blocks(whole[2000:], [1000, 3000]):
        grown.extend(block)
        filt.append(len(block))

    assert filt.ranges == ranges                    # the bounds are untouched
    reference = LocFilter(whole, **ranges)
    assert np.array_equal(filt.mask, reference.mask)


def test_a_hand_set_mask_excludes_new_rows_rather_than_guessing():
    grown = Localizations()
    grown.extend(_table(100))
    filt = LocFilter(grown)
    filt.set_mask("roi", np.ones(100, bool))
    grown.extend(_table(50, seed=1))
    filt.append(50)
    assert filt.mask[:100].all()
    assert not filt.mask[100:].any()


# ------------------------------------------------------------------ the state
@pytest.mark.parametrize("mode", ["hist", "precision"])
def test_an_appended_view_renders_what_the_whole_table_would(mode):
    whole = _table(20000)
    state = ViewState(Localizations(), RenderSettings(mode=mode), live=True,
                      extent=(0.0, 0.0, 10000.0, 8000.0))
    for block in _blocks(whole, [4000, 1, 5999, 10000]):
        state.append(block)
    state.filter.set("loc_precision_nm", None, 15.0)

    assert len(state.locs) == len(whole)
    for fov in (FieldOfView.from_range((0, 10000), (0, 8000), 20.0),
                FieldOfView.from_range((2000, 3000), (2000, 3000), 2.0)):
        reference = render_locs(whole, fov, state.settings,
                                select=LocFilter(
                                    whole, loc_precision_nm=(None, 15.0)).indices)
        assert np.allclose(state.render(fov).weight, reference.weight, atol=1e-6)


def test_appending_does_not_move_the_full_view_or_the_filter():
    state = ViewState(_table(3000), live=True,
                      extent=(0.0, 0.0, 10000.0, 8000.0))
    state.filter.set("z_nm", -150.0, 150.0)
    before_view = state.full_view()
    before_ranges = state.filter.ranges
    display, settings = state.display, state.settings

    state.append(_table(3000, seed=2, first_frame=3000))

    assert state.full_view() == before_view
    assert state.filter.ranges == before_ranges
    assert state.display is display and state.settings is settings
    assert len(state.filter) > 0


def test_the_grouped_table_goes_stale_rather_than_lagging_silently():
    state = ViewState(_table(2000), live=True)
    state.show_grouped(True)
    state.filter.set("loc_precision_nm", None, 12.0)   # the grouped filter
    grouped_before = state.locs
    assert not state.grouped_stale

    state.append(_table(2000, seed=3, first_frame=2000))
    assert state.grouped_stale
    assert state.locs is grouped_before                # the view does not jump

    rebuilt = state.group()
    assert not state.grouped_stale
    assert rebuilt.locs is not grouped_before
    assert rebuilt.filter.ranges["loc_precision_nm"] == (None, 12.0)


def test_appending_to_a_table_not_opened_for_it_is_refused():
    state = ViewState(_table(100))
    with pytest.raises(TypeError, match="live=True"):
        state.append(_table(10, seed=4))
