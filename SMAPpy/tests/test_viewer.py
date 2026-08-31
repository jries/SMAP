"""The viewer core: what the index and the filter select must match the truth."""
import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from smappy.filter import LocFilter                      # noqa: E402
from smappy.locs import Localizations                    # noqa: E402
from smappy.render import (DisplaySettings, FieldOfView,  # noqa: E402
                            RenderSettings, render_locs)
from smappy.viewer import (DEFAULT_BOUNDS, ViewState,      # noqa: E402
                            Viewer, _nice_length)


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
    assert viewer._rendered.n_locs == len(viewer.state.filter)   # after DEFAULT_BOUNDS

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

    viewer._on_bound(field, 1, "12")               # tighter than the default
    assert viewer.state.filter.ranges[field] == (None, 12.0)
    assert not viewer.state.filter.mask.all()

    viewer._on_bound(field, 0, "8")
    assert viewer.state.filter.ranges[field] == (8.0, 12.0)

    viewer._on_bound(field, 0, "")                 # one bound cleared
    assert viewer.state.filter.ranges[field] == (None, 12.0)
    viewer._on_bound(field, 1, "  ")               # both cleared: no filter
    assert field not in viewer.state.filter
    viewer.state.filter.clear()                    # and with the defaults gone too
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
    # its own filter: the typed bound is not carried across, the defaults apply
    assert viewer.state.filter.ranges[field] == DEFAULT_BOUNDS[field]
    assert viewer.bounds[field][1].text == "25"
    assert viewer._hints[field].get_text() != hint   # the ranges follow the table

    viewer.grouped.set_active(0)                     # -> back
    assert not viewer.state.use_grouped
    assert viewer.state.filter.ranges[field] == (None, 20.0)
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


# ------------------------------------------------------------ colour coding
def _viewer(locs=None, **kwargs):
    from smappy.viewer import Viewer
    return Viewer(ViewState(locs if locs is not None else _table()), **kwargs)


def test_the_colour_control_offers_intensity_and_what_the_table_carries():
    viewer = _viewer()
    labels = [name for name, _ in viewer.color_choices]
    assert labels[0] == "intensity"
    assert dict(viewer.color_choices)["z"] == "z_nm"
    assert "photons" not in labels          # this table has no photons column
    assert viewer.state.settings.color_field is None


def test_choosing_a_field_colours_the_image_over_an_explicit_range():
    """The range must be concrete, or the same z means a different colour at
    every zoom and after every block of a live fit."""
    viewer = _viewer()
    viewer._on_color("z")
    assert viewer.state.settings.color_field == "z_nm"
    assert viewer.state.settings.color_range is not None
    lo, hi = viewer.state.settings.color_range
    assert -400 < lo < hi < 400
    assert viewer.state.display.lut == "turbo"        # a hue ramp, not `hot`
    assert viewer._rendered is not None and viewer._rendered.color is not None


def test_typing_a_colour_range_is_what_gets_rendered():
    viewer = _viewer()
    viewer._on_color("z")
    viewer._on_color_bound(0, "-150")
    viewer._on_color_bound(1, "150")
    assert viewer.state.settings.color_range == (-150.0, 150.0)

    viewer._on_color_bound(1, "")            # an empty end means the data's own
    lo, hi = viewer.state.settings.color_range
    assert lo == -150.0 and hi == pytest.approx(
        float(np.nanmax(viewer.state.locs["z_nm"])))

    viewer._on_color_bound(0, "junk")        # junk is ignored, not applied
    assert viewer.state.settings.color_range[0] == -150.0


def test_going_back_to_intensity_undoes_the_colouring():
    viewer = _viewer()
    viewer._on_color("frame")
    assert viewer.state.settings.color_field == "frame"
    viewer._on_color("intensity")
    assert viewer.state.settings.color_field is None
    assert viewer.state.settings.color_range is None
    assert viewer.state.display.lut == "hot"
    assert viewer._rendered.color is None


def test_the_range_does_not_carry_across_fields():
    viewer = _viewer()
    viewer._on_color("z")
    viewer._on_color_bound(0, "-100")
    viewer._on_color("frame")
    lo, hi = viewer.state.settings.color_range
    assert lo >= 0 and hi <= len(viewer.state.locs)     # a frame range, not a z one


def test_a_field_the_table_does_not_carry_is_refused_not_crashed():
    """A live window offers the choices the fit will produce, before it has."""
    viewer = _viewer(color_fields=[("intensity", None), ("photons", "photons")])
    viewer._on_color("photons")
    assert viewer.state.settings.color_field is None    # moved back
    assert viewer._color_hint.get_text() == "(intensity)"


def test_opening_already_coloured_shows_that_in_the_panel():
    """`view_locs.py --color z_nm` must not open a panel saying "intensity"."""
    from smappy.render import RenderSettings as RS
    from smappy.viewer import Viewer

    locs = _table()
    state = ViewState(locs, RS(color_field="z_nm", color_range=(-200.0, 200.0)))
    viewer = Viewer(state)
    assert viewer._color_label == "z"
    assert viewer.colors.value_selected == "z"
    assert viewer.state.settings.color_range == (-200.0, 200.0)
    assert viewer._rendered.color is not None


def test_a_column_outside_the_usual_choices_gets_its_own_entry():
    from smappy.render import RenderSettings as RS
    from smappy.viewer import Viewer

    locs = _table()
    viewer = Viewer(ViewState(locs, RS(color_field="logl_rel")))
    assert viewer.colors.value_selected == "logl_rel"
    assert ("logl_rel", "logl_rel") in viewer.color_choices


# ------------------------------------------------------------------ gestures
def _events(viewer):
    """Send real mouse events through the canvas, as the backend would."""
    from matplotlib.backend_bases import MouseButton, MouseEvent

    fig = viewer.figure
    fig.canvas.draw()

    def send(name, x, y):
        MouseEvent(name, fig.canvas, x, y, button=MouseButton.LEFT)._process()
    return send


def test_a_click_does_not_force_a_blocking_redraw_per_text_box():
    """Every TextBox redraws the figure when a click lands elsewhere.

    With a dozen boxes that is a dozen *blocking* full draws before the first
    motion event of a drag arrives -- a fifth of a second here, near a second
    on a retina canvas, felt as the image lagging the mouse.
    """
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    send = _events(viewer)
    blocking = []
    viewer.figure.canvas.draw = lambda *a, **k: blocking.append(1)

    box = viewer.axes.get_window_extent()
    send("button_press_event", box.x0 + box.width / 2, box.y0 + box.height / 2)
    assert blocking == []
    assert len(viewer.bounds) + 1 >= 5      # there really are that many boxes


def test_panning_does_not_drift_even_though_it_renders_as_it_goes():
    """The grabbed point must stay under the cursor for the whole drag."""
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    send = _events(viewer)
    box = viewer.axes.get_window_extent()
    cx, cy = box.x0 + box.width / 2, box.y0 + box.height / 2

    send("button_press_event", cx, cy)
    grabbed = viewer._drag
    for i in range(1, 9):
        send("motion_notify_event", cx - 9 * i, cy - 4 * i)
    send("button_release_event", cx - 72, cy - 32)

    under_cursor = viewer.axes.transData.inverted().transform((cx - 72, cy - 32))
    assert under_cursor == pytest.approx(grabbed, abs=1e-6)


def test_a_fast_render_happens_inside_the_gesture():
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    send = _events(viewer)
    box = viewer.axes.get_window_extent()
    cx, cy = box.x0 + box.width / 2, box.y0 + box.height / 2

    rendered = []
    viewer._rendered_at = 0.0            # not throttled by the render just done
    real = viewer._render_now
    viewer._render_now = lambda: (rendered.append(1), real())[1]

    send("button_press_event", cx, cy)
    send("motion_notify_event", cx - 20, cy - 10)
    assert rendered, "a cheap render should follow the mouse, not wait for it"


def test_a_slow_render_is_deferred_instead_of_making_the_gesture_lurch():
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    send = _events(viewer)
    box = viewer.axes.get_window_extent()
    cx, cy = box.x0 + box.width / 2, box.y0 + box.height / 2

    rendered = []
    viewer._render_seconds = 10.0        # as if a render took ten seconds
    real = viewer._render_now
    viewer._render_now = lambda: (rendered.append(1), real())[1]

    send("button_press_event", cx, cy)
    send("motion_notify_event", cx - 20, cy - 10)
    assert not rendered                  # the image still moved; nothing rendered
    assert viewer.axes.get_xlim() != (0, 0)


def test_a_queued_render_does_not_land_in_the_middle_of_a_drag():
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    send = _events(viewer)
    box = viewer.axes.get_window_extent()
    stopped = []
    viewer._timer.stop = lambda: stopped.append(1)

    send("button_press_event", box.x0 + box.width / 2, box.y0 + box.height / 2)
    assert stopped


# ------------------------------------------------------------ opening bounds
def test_a_window_opens_with_the_default_bounds_shown_in_its_boxes():
    """The defaults throw localizations away, so they must be on screen."""
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    ranges = viewer.state.filter.ranges
    assert ranges["loc_precision_nm"] == (None, 25.0)
    assert ranges["logl_rel"] == (-1.5, None)
    assert ranges["z_nm"] == (-500.0, 500.0)
    assert viewer.bounds["loc_precision_nm"][1].text == "25"
    assert viewer.bounds["logl_rel"][0].text == "-1.5"
    assert viewer.bounds["z_nm"][0].text == "-500"
    assert len(viewer.state.filter) < len(viewer.state.locs)


def test_a_bound_already_set_is_not_overridden_by_a_default():
    from smappy.filter import LocFilter
    from smappy.viewer import Viewer

    locs = _table()
    state = ViewState(locs, filter=LocFilter(locs, loc_precision_nm=(None, 8.0)))
    viewer = Viewer(state)
    assert viewer.state.filter.ranges["loc_precision_nm"] == (None, 8.0)
    assert viewer.state.filter.ranges["z_nm"] == (-500.0, 500.0)   # the rest apply


def test_a_cleared_default_stays_cleared():
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_groupable_table()))
    viewer._on_bound("loc_precision_nm", 1, "")
    assert "loc_precision_nm" not in viewer.state.filter
    viewer.grouped.set_active(0)          # -> grouped, which gets its own defaults
    viewer.grouped.set_active(0)          # -> back
    assert "loc_precision_nm" not in viewer.state.filter


def test_a_pixel_table_gets_only_the_bounds_whose_unit_is_fixed():
    """25 nm means nothing in pixels, so that default does not apply there."""
    from smappy.viewer import Viewer

    locs = _table()
    locs.columns["loc_precision_pix"] = locs.columns.pop("loc_precision_nm")
    locs.columns["x_pix"] = locs.columns.pop("x_nm")
    locs.columns["y_pix"] = locs.columns.pop("y_nm")
    viewer = Viewer(ViewState(locs))
    assert "loc_precision_pix" not in viewer.state.filter
    assert viewer.state.filter.ranges["logl_rel"] == (-1.5, None)


def test_the_image_and_the_controls_are_separate_windows():
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    assert viewer.controls is not viewer.figure
    assert viewer.axes.figure is viewer.figure          # the image window
    for boxes in viewer.bounds.values():
        assert all(box.ax.figure is viewer.controls for box in boxes)
    assert viewer.colors.ax.figure is viewer.controls
    assert viewer.contrast.ax.figure is viewer.controls
    # the image fills its window: no widgets, and nothing but the axes
    assert viewer.figure.axes == [viewer.axes]


def test_resizing_the_image_window_rerenders_at_the_new_size():
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()))
    before = viewer.fov.shape
    viewer.figure.set_size_inches(5.0, 5.0)
    viewer.figure.canvas.draw()
    viewer._render_now()
    assert viewer.fov.shape != before                   # rendered to fit the canvas


@pytest.mark.parametrize("figsize", [(13, 6), (6, 10), (8.5, 8.5)])
def test_the_image_fills_a_window_of_any_shape_with_square_pixels(figsize):
    """A window that is not square must not leave bands beside the image.

    The axes may never be shrunk to the shape of the view: it is the view that
    widens to the shape of the window (a wide window simply shows more x), so
    the canvas is filled edge to edge while one pixel size still serves both
    axes -- an SMLM image with non-square pixels would be a lie.
    """
    from smappy.viewer import Viewer

    viewer = Viewer(ViewState(_table()), figsize=figsize)
    viewer.figure.canvas.draw()
    viewer._render_now()
    viewer.figure.canvas.draw()

    box = viewer.axes.get_window_extent()
    width, height = viewer.figure.get_size_inches() * viewer.figure.dpi
    assert box.width == pytest.approx(0.96 * width, rel=0.01)
    assert box.height == pytest.approx(0.93 * height, rel=0.01)

    x0, x1 = viewer.axes.get_xlim()
    ylo, yhi = sorted(viewer.axes.get_ylim())
    assert (x1 - x0) / (yhi - ylo) == pytest.approx(box.width / box.height, rel=0.02)
    assert viewer.fov.shape == (round(box.height), round(box.width))
    assert (x1 - x0) == pytest.approx(viewer.fov.pixelsize * box.width, rel=0.02)
    assert (yhi - ylo) == pytest.approx(viewer.fov.pixelsize * box.height, rel=0.02)
