"""Cached range filters, and rendering a table through them."""
import numpy as np
import pytest

from smappy.filter import LocFilter, quantile_range
from smappy.locs import Localizations
from smappy.render import (DisplaySettings, FieldOfView, RenderSettings,
                            render_locs)


def _table(n=1000, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.uniform(-400, 400, n).astype(np.float32)
    z[::50] = np.nan                       # localizations without a z
    return Localizations({
        "x_nm": rng.uniform(0, 10000, n).astype(np.float32),
        "y_nm": rng.uniform(0, 10000, n).astype(np.float32),
        "z_nm": z,
        "loc_precision_nm": rng.gamma(4, 3, n).astype(np.float32),
        "logl_rel": rng.normal(-1, 0.5, n).astype(np.float32),
        "frame": np.arange(n, dtype=np.int64),
    }, {"units": "nm"})


def test_ranges_intersect():
    locs = _table()
    f = LocFilter(locs, loc_precision_nm=(None, 15.0), logl_rel=(-2.0, None))
    expected = (locs["loc_precision_nm"] <= 15.0) & (locs["logl_rel"] >= -2.0)
    assert np.array_equal(f.mask, expected)
    assert len(f) == int(expected.sum()) and 0 < len(f) < len(locs)
    assert np.array_equal(f.indices, np.flatnonzero(expected))


def test_nan_is_excluded_by_any_range():
    """A localization with no z must not appear in a z slab."""
    locs = _table()
    missing = np.isnan(locs["z_nm"])
    assert missing.any()
    f = LocFilter(locs, z_nm=(-1e9, 1e9))
    assert not f.mask[missing].any()
    # ...but an unbounded filter keeps everything, NaN included
    assert LocFilter(locs).mask.all()


def test_changing_one_filter_recomputes_only_that_one():
    locs = _table()
    f = LocFilter(locs, loc_precision_nm=(0.0, 15.0), logl_rel=(-2.0, 0.0))
    kept = f.mask.copy()
    other = f._masks["logl_rel"]           # the array a GUI must not have to redo

    f.set("loc_precision_nm", 0.0, 10.0)
    assert f._masks["logl_rel"] is other   # untouched, not merely equal
    assert f.mask.sum() < kept.sum()       # and the result did change

    f.remove("loc_precision_nm")
    assert np.array_equal(f.mask, other)
    assert f.clear().mask.all()


def test_arbitrary_masks_join_the_intersection():
    locs = _table()
    inside = locs["x_nm"] < 5000
    f = LocFilter(locs, loc_precision_nm=(None, 15.0)).set_mask("view", inside)
    assert np.array_equal(f.mask, inside & (locs["loc_precision_nm"] <= 15.0))
    assert f.counts()["view"] == int(inside.sum())
    with pytest.raises(ValueError):
        f.set_mask("bad", np.ones(3, bool))
    with pytest.raises(KeyError):
        f.set("not_a_column", 0, 1)


def test_quantile_range_ignores_tails_and_nan():
    locs = _table()
    lo, hi = quantile_range(locs, "z_nm")
    assert np.isfinite([lo, hi]).all() and -400 < lo < hi < 400


def test_filter_and_render_agree_with_rendering_the_filtered_table():
    locs = _table()
    fov = FieldOfView.from_range((0, 10000), (0, 10000), 50.0)
    f = LocFilter(locs, loc_precision_nm=(None, 12.0))
    settings = RenderSettings(mode="precision")

    through_filter = render_locs(locs, fov, settings, select=f)
    on_subset = render_locs(locs[f.mask], fov, settings)
    assert np.array_equal(through_filter.weight, on_subset.weight)
    assert through_filter.n_locs == len(f)
    # and it is genuinely a subset of the unfiltered image
    assert through_filter.weight.sum() < render_locs(locs, fov, settings).weight.sum()


def test_field_colouring_uses_the_display_lut():
    locs = _table()
    fov = FieldOfView.from_range((0, 10000), (0, 10000), 50.0)
    f = LocFilter(locs, z_nm=(-400, 400))          # drops the NaN z values
    settings = RenderSettings(mode="hist", color_field="z_nm", color_range=(-400, 400))

    img = render_locs(locs, fov, settings, DisplaySettings(lut="jet"), select=f)
    assert img.is_colored and img.n_locs == len(f)
    rgb = DisplaySettings(lut="jet").apply(img)
    assert rgb.shape == (fov.ny, fov.nx, 3) and rgb.max() <= 1.0
    # an intensity render of the same selection carries the same weight plane
    plain = render_locs(locs, fov, RenderSettings(mode="hist"), select=f)
    assert np.array_equal(plain.weight, img.weight)
