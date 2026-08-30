"""The viewer core: what the index and the filter select must match the truth."""
import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from smapfit.filter import LocFilter                      # noqa: E402
from smapfit.locs import Localizations                    # noqa: E402
from smapfit.render import (DisplaySettings, FieldOfView,  # noqa: E402
                            RenderSettings, render_locs)
from smapfit.viewer import ViewState, Viewer, _nice_length  # noqa: E402


def _table(n=20000, seed=0):
    rng = np.random.default_rng(seed)
    return Localizations({
        "x_nm": rng.uniform(0, 10000, n).astype(np.float32),
        "y_nm": rng.uniform(0, 8000, n).astype(np.float32),
        "z_nm": rng.uniform(-400, 400, n).astype(np.float32),
        "loc_precision_nm": (rng.gamma(4, 3, n) + 5).astype(np.float32),
        "logl_rel": rng.normal(-1, 0.5, n).astype(np.float32),
        "frame": np.arange(n, dtype=np.int64),
    }, {"units": "nm"})


@pytest.mark.parametrize("mode", ["hist", "precision"])
def test_the_index_does_not_change_the_image(mode):
    """Selecting through the spatial index must match rendering everything.

    This is the whole safety property of the index: it may hand back extra
    candidates, but it may never lose one whose blob reaches the view.
    """
    locs = _table()
    state = ViewState(locs, RenderSettings(mode=mode))
    state.filter.set("loc_precision_nm", None, 15.0)
    for fov in (FieldOfView.from_range((0, 10000), (0, 8000), 20.0),
                FieldOfView.from_range((2000, 3000), (2000, 3000), 2.0),
                FieldOfView.from_range((-5000, -4000), (0, 1000), 2.0)):
        through_index = state.render(fov)
        everything = render_locs(locs, fov, state.settings,
                                 select=state.filter.indices)
        assert np.allclose(through_index.weight, everything.weight, atol=1e-6)


def test_a_filter_change_is_visible_without_rebuilding_the_index():
    locs = _table()
    state = ViewState(locs)
    index = state.index
    fov = FieldOfView.from_range((0, 10000), (0, 8000), 20.0)
    before = state.render(fov).n_locs
    state.filter.set("loc_precision_nm", None, 12.0)
    assert state.render(fov).n_locs < before
    assert state.index is index          # positions did not change, so it stands


def test_viewer_zoom_pan_and_reset():
    viewer = Viewer(ViewState(_table()))
    full = viewer.fov
    assert viewer._rendered.n_locs == len(viewer.state.locs)

    viewer._zoom(0.25, cx=5000.0, cy=4000.0)
    viewer._render_now()
    assert viewer.fov.pixelsize == pytest.approx(full.pixelsize * 0.25, rel=0.02)
    assert viewer._rendered.n_locs < full.pixelsize * 0  + len(viewer.state.locs)
    assert viewer.fov.shape == full.shape          # the canvas is unchanged

    x0 = viewer.fov.x0
    viewer.axes.set_xlim(x0 + 500, viewer.fov.x1 + 500)
    viewer._render_now()
    assert viewer.fov.x0 == pytest.approx(x0 + 500, abs=viewer.fov.pixelsize)

    viewer.reset()
    assert viewer.fov.pixelsize == pytest.approx(full.pixelsize, rel=0.02)


def test_typed_bounds_filter_and_an_empty_box_means_no_bound():
    viewer = Viewer(ViewState(_table()))
    field = "loc_precision_nm"

    viewer._on_bound(field, 1, "12")
    assert viewer.state.filter.ranges[field] == (None, 12.0)
    assert not viewer.state.filter.mask.all()

    viewer._on_bound(field, 0, "8")
    assert viewer.state.filter.ranges[field] == (8.0, 12.0)

    viewer._on_bound(field, 0, "")                 # one bound cleared
    assert viewer.state.filter.ranges[field] == (None, 12.0)
    viewer._on_bound(field, 1, "  ")               # both cleared: no filter
    assert field not in viewer.state.filter
    assert viewer.state.filter.mask.all()


def test_junk_in_a_box_changes_nothing_and_is_put_back():
    viewer = Viewer(ViewState(_table()))
    viewer._on_bound("loc_precision_nm", 1, "12")
    viewer._on_bound("loc_precision_nm", 1, "not a number")
    assert viewer.state.filter.ranges["loc_precision_nm"] == (None, 12.0)
    assert viewer.bounds["loc_precision_nm"][1].text == "12"


def test_contrast_and_gamma_redisplay_without_rendering():
    viewer = Viewer(ViewState(_table()))
    rendered = viewer._rendered
    before = viewer.image.get_array().copy()
    viewer._on_contrast(2.0)                       # brighter
    viewer._on_gamma(0.5)
    assert viewer._rendered is rendered            # the planes were reused
    assert not np.array_equal(before, viewer.image.get_array())


def test_scalebar_length_is_a_round_number():
    for span in (10000.0, 1234.0, 0.9, 33.0):
        length = _nice_length(span)
        assert length <= span / 5.0 < 2 * length   # a fifth of the view, roughly
        mantissa = length / 10.0 ** np.floor(np.log10(length))
        assert mantissa == pytest.approx(round(mantissa)) and mantissa in (1, 2, 5)


def _groupable_table(n=3000, seed=5):
    """A table where each emitter is on for a few consecutive frames."""
    rng = np.random.default_rng(seed)
    emitters = rng.uniform(0, 5000, (n // 3, 2))
    x, y, frame = [], [], []
    for i, (px, py) in enumerate(emitters):
        for k in range(3):
            x.append(px + rng.normal(0, 5)); y.append(py + rng.normal(0, 5))
            frame.append(i % 500 + k)
    n = len(x)
    return Localizations({
        "x_nm": np.array(x, np.float32), "y_nm": np.array(y, np.float32),
        "frame": np.array(frame, np.int64),
        "loc_precision_nm": (rng.gamma(4, 3, n) + 5).astype(np.float32),
        "photons": rng.gamma(3, 500, n).astype(np.float32),
        "logl_rel": rng.normal(-1, 0.5, n).astype(np.float32),
    }, {"units": "nm"})


def test_grouped_and_ungrouped_keep_separate_filters():
    state = ViewState(_groupable_table())
    state.filter.set("loc_precision_nm", None, 20.0)

    state.show_grouped(True)
    assert state.has_grouped
    assert len(state.locs) < len(state.sets["ungrouped"].locs)
    assert state.filter.ranges == {}                 # its own, untouched filter
    state.filter.set("loc_precision_nm", None, 8.0)

    state.show_grouped(False)
    assert state.filter.ranges == {"loc_precision_nm": (None, 20.0)}
    state.show_grouped(True)
    assert state.filter.ranges == {"loc_precision_nm": (None, 8.0)}


def test_switching_does_not_regroup():
    state = ViewState(_groupable_table())
    state.show_grouped(True)
    grouped = state.sets["grouped"]
    state.show_grouped(False)
    state.show_grouped(True)
    assert state.sets["grouped"] is grouped          # built once, then reused


def test_each_set_has_its_own_index_and_renders():
    state = ViewState(_groupable_table())
    fov = FieldOfView.from_range((0, 5000), (0, 5000), 20.0)
    ungrouped = state.render(fov).n_locs
    state.show_grouped(True)
    assert state.index is not state.sets["ungrouped"].index
    assert state.render(fov).n_locs < ungrouped


def test_the_viewer_switch_restores_the_boxes_silently():
    viewer = Viewer(ViewState(_groupable_table()))
    field = "loc_precision_nm"
    viewer._on_bound(field, 1, "20")
    hint = viewer._hints[field].get_text()

    viewer.grouped.set_active(0)                     # -> grouped
    assert viewer.state.use_grouped
    assert viewer.state.filter.ranges == {}          # not carried across
    assert viewer.bounds[field][1].text == ""
    assert viewer._hints[field].get_text() != hint   # the ranges follow the table

    viewer.grouped.set_active(0)                     # -> back
    assert not viewer.state.use_grouped
    assert viewer.state.filter.ranges == {field: (None, 20.0)}
    assert viewer.bounds[field][1].text == "20"
    assert viewer._hints[field].get_text() == hint


def test_the_additive_switch_only_redisplays():
    locs = _table()
    state = ViewState(locs, RenderSettings(mode="hist", color_field="z_nm",
                                           color_range=(-400, 400)),
                      DisplaySettings(lut="turbo"))
    viewer = Viewer(state)
    rendered = viewer._rendered
    assert state.display.color_mode == "hue"
    before = viewer.image.get_array().copy()

    viewer.switches.set_active(1)                    # -> additive
    assert state.display.color_mode == "sum"
    assert viewer._rendered is rendered              # the planes were reused
    assert not np.array_equal(before, viewer.image.get_array())

    viewer.switches.set_active(1)
    assert state.display.color_mode == "hue"
    assert np.allclose(before, viewer.image.get_array())
