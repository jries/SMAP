"""An interactive view of a localization table: zoom, pan, filter.

Two layers.  :class:`ViewState` is the part that has nothing to do with a GUI:
it owns the table, the spatial index, the filter and the settings, and turns a
field of view into an RGB image.  :class:`Viewer` is a matplotlib front end on
top of it.  Anything else -- Qt, a notebook widget -- reuses the first and
replaces the second.

Three decisions make this usable on millions of localizations:

* **Render at display resolution.**  The pixel size follows from the window, not
  from the zoom, so the cost of a render is bounded by the canvas: zoomed out
  the blobs collapse onto the `mingausspix` floor, zoomed in there are fewer
  localizations on screen.
* **Gestures are image-space; re-rendering is deferred.**  A drag or a scroll
  only moves the axes limits, and matplotlib rescales the image it already has;
  a timer re-renders once the gesture has been still for a moment.  So the
  interaction stays smooth even where a render takes 130 ms.
* **The spatial index answers what is on screen**, so a zoomed view never walks
  the whole table.

The filter controls cover the fields that are actually used -- localization
precision, z (or the PSF size for a free fit), relative log-likelihood, frame --
rather than every column, and each is a typed minimum and maximum rather than a
slider: the numbers wanted are exact ones (a 20 nm precision cut, a +-150 nm z
slab), and the useful ranges span decades.  An empty box means no bound.  A
window opens with the few bounds in `DEFAULT_BOUNDS` -- the ones that are right
almost every time -- and each of them is written into its box, so what has been
filtered out is on screen rather than hidden in a default.  The data range is
printed beside each row.

The image and the controls are separate windows: the image can then be as large
as the screen allows, and because the rendered pixel size follows the canvas,
making it larger renders a finer image rather than scaling one up.
"""

from __future__ import annotations

import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .filter import LocFilter, quantile_range
from .group import GroupSettings, group
from .locs import Localizations
from .render import (DisplaySettings, FieldOfView, RenderSettings, RenderedImage,
                     positions, render_locs)
from .spatial import GrowingIndex, SpatialIndex

# The fields worth a filter control, best alternative first.  Everything else
# stays available through `LocFilter` without cluttering the window.
FILTER_FIELDS: Sequence[Tuple[str, ...]] = (
    ("loc_precision_nm", "loc_precision_pix"),
    ("z_nm",),
    ("sigma_nm", "sigma_pix"),
    ("logl_rel",),
    ("frame",),
)

# What an image can be coloured by, beside plain intensity: a short label and
# the column it stands for, best alternative first.  Colour is baked into the
# accumulation (each localization contributes its own colour), so changing this
# is a re-render, unlike contrast or gamma.
COLOR_FIELDS: Sequence[Tuple[str, Tuple[str, ...]]] = (
    ("z", ("z_nm",)),
    ("frame", ("frame",)),
    ("precision", ("loc_precision_nm", "loc_precision_pix")),
    ("photons", ("photons",)),
)

# The bounds a window opens with.  Every one of them throws localizations
# away, so each is shown in its box: what is filtered out is on screen, not
# hidden in a default.  Only fields whose unit is fixed get one -- a 25 nm
# precision cut means nothing in pixels -- and a bound already set (a file
# opened with one, a table switched back to) is never overridden.
DEFAULT_BOUNDS: Dict[str, Tuple[Optional[float], Optional[float]]] = {
    "loc_precision_nm": (None, 25.0),   # nm; a looser cut is rarely wanted
    "logl_rel": (-1.5, None),           # drops the fits that did not converge
    "z_nm": (-500.0, 500.0),            # the depth a spline calibration covers
}

INTENSITY_LUT = "hot"    # an intensity ramp: brightness is the only dimension
FIELD_LUT = "turbo"      # a full-hue ramp, which is what a coded field needs


class LocSet:
    """One table with the things that are derived from it and nothing else.

    Grouped and ungrouped localizations are different tables, so each needs its
    own filter and its own spatial index.  Keeping both alive is the point:
    linking costs a third of a second and the filter caches would have to be
    rebuilt on every switch otherwise.
    """

    def __init__(self, locs: Localizations, filter: Optional[LocFilter] = None,
                 extent: Optional[Sequence[float]] = None,
                 growable: bool = False):
        self.locs = locs
        self.filter = filter if filter is not None else LocFilter(locs)
        # a live view is opened before the first block has been fitted, so an
        # empty table is a normal starting point, not an error
        empty = np.empty(0, np.float32)
        x, y = positions(locs) if locs.columns else (empty, empty)
        self.index = (GrowingIndex(x, y, extent=extent) if growable
                      else SpatialIndex(x, y, extent=extent))

        # for the query margin: how far the widest blob reaches beyond the view
        self.median_precision = 0.0
        self._precision_field: Optional[str] = None
        for name in ("loc_precision_nm", "loc_precision_pix"):
            if name in locs:
                self._precision_field = name
                break
        self._precision_n = 0
        self._update_precision()

    def _update_precision(self, force: bool = True) -> None:
        """Recompute the median localization precision.

        It only sets the cap on the rendering sigma and it hardly moves when a
        few percent more localizations arrive, so while data is streaming in it
        is recomputed once the table has grown by a fifth rather than on every
        block -- a median over millions of values is not free.
        """
        if self._precision_field is None:
            return
        if not force and len(self.locs) < self._precision_n * 1.2:
            return
        values = np.asarray(self.locs[self._precision_field], np.float64)
        values = values[np.isfinite(values)]
        if values.size:
            self.median_precision = float(np.median(values))
        self._precision_n = len(self.locs)

    def append(self, new: Localizations) -> int:
        """Add localizations to the end of the table.  Returns how many.

        The table, the filter and the index are all extended rather than
        rebuilt, and nothing that was set -- a filter bound, the grid the index
        is on -- is re-derived from the larger table.
        """
        n = len(new)
        if n == 0:
            return 0
        if not isinstance(self.index, GrowingIndex):
            raise TypeError("this table was not opened for appending; "
                            "build the ViewState with live=True")
        x, y = positions(new)
        self.locs.extend(new)
        self.filter.append(n)
        self.index.append(x, y)
        self._update_precision(force=False)
        return n

    def __len__(self) -> int:
        return len(self.locs)


class ViewState:
    """Everything needed to render a view, without a window around it."""

    def __init__(self, locs: Localizations,
                 settings: Optional[RenderSettings] = None,
                 display: Optional[DisplaySettings] = None,
                 filter: Optional[LocFilter] = None,
                 grouped: Optional[Localizations] = None,
                 live: bool = False,
                 extent: Optional[Sequence[float]] = None):
        """``live`` prepares the table for `append`; ``extent`` fixes the area
        the index covers, so the full view does not move as data arrives."""
        self.settings = settings or RenderSettings()
        self.display = display or DisplaySettings()
        self.n_threads = 0
        self.sets: Dict[str, LocSet] = {
            "ungrouped": LocSet(locs, filter, extent=extent, growable=live)}
        if grouped is not None:
            self.sets["grouped"] = LocSet(grouped)
        self.use_grouped = False
        self.grouped_stale = False

    # The active table.  Everything downstream goes through these, so switching
    # is a single assignment and no cached mask is thrown away.
    @property
    def current(self) -> LocSet:
        return self.sets["grouped" if self.use_grouped else "ungrouped"]

    @property
    def locs(self) -> Localizations:
        return self.current.locs

    @property
    def filter(self) -> LocFilter:
        return self.current.filter

    @property
    def index(self) -> SpatialIndex:
        return self.current.index

    @property
    def has_grouped(self) -> bool:
        return "grouped" in self.sets

    def group(self, settings: Optional[GroupSettings] = None) -> LocSet:
        """Build the grouped table, once.  This is the slow part.

        Grouping links localizations across frames and has no incremental form,
        so an append cannot extend it: the table is marked stale instead and
        rebuilt the next time it is asked for.
        """
        if "grouped" not in self.sets or self.grouped_stale:
            grouped, _ = group(self.sets["ungrouped"].locs, settings)
            keep = self.sets.get("grouped")
            self.sets["grouped"] = LocSet(grouped)
            if keep is not None:      # a rebuild: the bounds the user set stand
                for field, (lo, hi) in keep.filter.ranges.items():
                    if field in grouped:
                        self.sets["grouped"].filter.set(field, lo, hi)
            self.grouped_stale = False
        return self.sets["grouped"]

    def append(self, new: Localizations) -> int:
        """Add localizations to the ungrouped table.  Returns how many.

        Called from whichever thread owns the view -- the GUI one -- so that
        nothing here has to be thread-safe; a fitter running in the background
        hands blocks over through a queue.
        """
        n = self.sets["ungrouped"].append(new)
        if n and "grouped" in self.sets:
            self.grouped_stale = True
        return n

    def show_grouped(self, on: bool,
                     settings: Optional[GroupSettings] = None) -> None:
        if on:
            self.group(settings)
        self.use_grouped = bool(on)

    def full_view(self, margin_fraction: float = 0.01) -> Tuple[tuple, tuple]:
        """The coordinate ranges covering every localization."""
        x0, y0, x1, y1 = self.index.bounds
        mx, my = (x1 - x0) * margin_fraction, (y1 - y0) * margin_fraction
        return (x0 - mx, x1 + mx), (y0 - my, y1 + my)

    def max_sigma(self, pixelsize: float) -> float:
        """The largest rendering sigma at this zoom, in data units."""
        if self.settings.mode == "hist":
            return 0.0
        if self.settings.mode == "gauss":
            return self.settings.sigma
        s = self.settings.sigma_settings
        floor = max(s.min_sigma, s.min_sigma_pixels * pixelsize)
        return max(floor, s.max_factor * self.current.median_precision * s.factor)

    def select(self, fov: FieldOfView) -> np.ndarray:
        """The localizations that can contribute to this view, after filtering."""
        margin = self.settings.roi_sigma * self.max_sigma(fov.pixelsize)
        candidates = self.index.query(fov.x0, fov.x1, fov.y0, fov.y1, margin)
        if candidates.size == self.index.n_localizations:
            return self.filter.indices          # nothing to narrow: skip the gather
        return candidates[self.filter.mask[candidates]]

    def render(self, fov: FieldOfView) -> RenderedImage:
        return render_locs(self.locs, fov, self.settings, self.display,
                           select=self.select(fov), n_threads=self.n_threads)

    def image(self, fov: FieldOfView) -> Tuple[np.ndarray, RenderedImage]:
        rendered = self.render(fov)
        return self.display.apply(rendered), rendered


def _nice_length(span: float) -> float:
    """A round number about a fifth of the span, for the scale bar."""
    raw = span / 5.0
    power = 10.0 ** np.floor(np.log10(max(raw, 1e-12)))
    return max((step * power for step in (1.0, 2.0, 5.0) if step * power <= raw),
               default=power)


def _window_title(figure, title: str) -> None:
    """Name a window, where the backend has windows to name (Agg has not)."""
    manager = getattr(figure.canvas, "manager", None)
    if manager is not None and hasattr(manager, "set_window_title"):
        manager.set_window_title(title)


def _patched_text_box(TextBox):
    """A `TextBox` that does not force a blocking redraw of the whole figure.

    Every TextBox redraws the figure when a click lands anywhere else, to take
    its cursor away -- with `canvas.draw()`, which blocks.  A window with a
    dozen boxes therefore pays a dozen full draws before the first motion event
    of a drag is even delivered: a fifth of a second here, close to a second on
    a retina canvas, felt as the image not moving when the mouse does.

    `draw_idle` collapses all of them into the one draw that was going to
    happen anyway.  It is swapped in around the call rather than reimplementing
    `stop_typing`, whose body reaches into private state that changes between
    matplotlib versions.
    """

    class _TextBox(TextBox):
        def stop_typing(self):
            canvas = self.ax.figure.canvas
            had = canvas.__dict__.get("draw")
            canvas.draw = canvas.draw_idle
            try:
                super().stop_typing()
            finally:
                if had is None:
                    canvas.__dict__.pop("draw", None)
                else:
                    canvas.draw = had

    return _TextBox


class Viewer:
    """Two matplotlib windows over a :class:`ViewState`: an image and controls.

    The image has a window to itself so it can be resized to whatever the
    screen allows -- and since the rendered pixel size follows the canvas, a
    bigger window is a *finer* image rather than a scaled-up one.  The controls
    sit in a second, small window that keeps its size, so they neither eat into
    the image nor stretch across it.  Closing the image closes both; closing
    the controls leaves the image alone.

    Scroll or pinch to zoom about the cursor, ``+``/``-`` to zoom about the
    centre, drag to pan, ``r`` to reset.  Each filter field takes a minimum and
    a maximum, applied on Enter; an empty box means no bound.  The two sliders
    set the contrast and gamma, which re-display without re-rendering.
    "grouped" switches to the grouped table, building it on first use;
    "additive" switches the field-colour composite to SMAP's.
    """

    REDRAW_DELAY_MS = 150
    # A render this fast can be run inside the gesture itself, so panning and
    # zooming fill in what they expose instead of dragging a stale image around
    # and only catching up once the mouse stops.  Above it, moving the image is
    # all that happens until the gesture settles: a slow render inside the
    # gesture would make it lurch, which is worse than a blank edge.
    LIVE_RENDER_SECONDS = 0.12

    def __init__(self, state: ViewState, figsize=(8.5, 8.5),
                 filter_fields: Optional[Sequence[str]] = None,
                 color_fields: Optional[Sequence[Tuple[str, Optional[str]]]] = None,
                 group_settings: Optional[GroupSettings] = None,
                 control_width: float = 4.6):
        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, RadioButtons, Slider
        from matplotlib.widgets import TextBox as _MplTextBox

        TextBox = _patched_text_box(_MplTextBox)

        self.state = state
        self.group_settings = group_settings
        self.fov: Optional[FieldOfView] = None
        self._rendered: Optional[RenderedImage] = None
        self.status = ""            # appended to the title; the live view uses it
        self._closed = False

        fields = list(filter_fields) if filter_fields is not None \
            else self._available_fields()
        self.color_choices: List[Tuple[str, Optional[str]]] = (
            list(color_fields) if color_fields is not None
            else self._available_color_fields())

        # The image gets the whole of its own window, so it can be made as
        # large as the screen allows and the render follows: the rendered pixel
        # size comes from the canvas, so a bigger window is a finer image
        # rather than a scaled-up one.
        self.figure = plt.figure(figsize=figsize)
        _window_title(self.figure, "smappy")
        self.axes = self.figure.add_axes([0.02, 0.02, 0.96, 0.93])
        self.axes.set_facecolor("black")
        # Square pixels, but never at the cost of the window: with
        # ``adjustable="datalim"`` matplotlib widens the view to fill the axes
        # instead of shrinking the axes to the view's shape.  So a window of
        # any proportion is filled edge to edge -- a wide one simply shows more
        # x -- and `FieldOfView.fit`, which sizes a render from the canvas and
        # keeps one pixel size for both axes, agrees with it by construction.
        self.axes.set_aspect("equal", adjustable="datalim")

        # The controls get a second, small window.  Rows are laid out in inches
        # and the window is sized to fit them, so the spacing stays the same
        # whatever the number of filter fields or colour choices.
        row, box_h, slider_h, gap, pad = 0.34, 0.24, 0.16, 0.12, 0.16
        block = max(0.95, 0.20 * len(self.color_choices) + 0.30)
        height = (pad + block + 0.34 + 2 * (slider_h + 0.16) + gap
                  + (len(fields) + 1) * row + pad)
        self.controls = plt.figure(figsize=(control_width, height))
        _window_title(self.controls, "smappy controls")

        def place(x, y_inches, w, h):
            """An axes in the control window, positioned in inches from below."""
            return self.controls.add_axes([x, y_inches / height, w, h / height])

        def hint(x, y_inches):
            return self.controls.text(x, (y_inches + 0.07) / height, "",
                                      fontsize=7, color="0.45")

        y = pad
        # "grouped" rebuilds the table (once); "additive" only changes how the
        # colour planes are composited, so it re-displays without re-rendering
        ax = place(0.52, y, 0.44, block)
        ax.set_frame_on(False)
        self.switches = CheckButtons(
            ax, ["grouped", "additive"],
            [state.use_grouped, state.display.color_mode == "sum"])
        for text in self.switches.labels:
            text.set_fontsize(8)
        self.switches.on_clicked(self._on_switch)
        self.grouped = self.switches   # the switch panel, kept under both names

        # What to colour by.  Mutually exclusive, so radio buttons rather than
        # another check box: an image is coloured by one thing or by nothing.
        # The panel starts on whatever the settings passed in already say, so
        # opening with a colour field selected does not show "intensity".
        active = self._color_choice_for(self.state.settings.color_field)
        ax = place(0.08, y, 0.38, block)
        ax.set_frame_on(False)
        ax.set_title("colour by", fontsize=8, color="0.3")
        self.colors = RadioButtons(ax, [name for name, _ in self.color_choices],
                                   active=active)
        for text in self.colors.labels:
            text.set_fontsize(8)
        self.colors.on_clicked(self._on_color)
        self._color_label = self.color_choices[active][0]
        y += block + 0.34       # room for the "colour by" title above the panel

        # SMAP's dynamic contrast: saturate 10^-p of the pixels.  The useful
        # range spans decades, which is why the slider moves the exponent.
        self.gamma = Slider(place(0.30, y, 0.50, slider_h), "gamma", 0.2, 2.0,
                            valinit=state.display.gamma)
        self.gamma.on_changed(self._on_gamma)
        y += slider_h + 0.16
        self.contrast = Slider(place(0.30, y, 0.50, slider_h), "contrast (p)",
                               1.0, 5.0, valinit=state.display.contrast)
        self.contrast.on_changed(self._on_contrast)
        for slider in (self.contrast, self.gamma):
            slider.label.set_fontsize(8)
            slider.valtext.set_fontsize(8)
        y += slider_h + 0.16 + gap

        # Colouring by a field needs a range, and it must be an explicit one:
        # taken from the data it would be re-derived on every render, so the
        # same z would mean a different colour at every zoom and after every
        # block of a live fit.  The row reads as the filter rows do -- a
        # minimum and a maximum, with the field and its data range beside them.
        self._color_typed: List[Optional[float]] = [None, None]
        low = TextBox(place(0.40, y, 0.17, box_h), "colour", initial="",
                      textalignment="right")
        high = TextBox(place(0.60, y, 0.17, box_h), "to", initial="",
                       textalignment="right")
        for box, which in ((low, 0), (high, 1)):
            box.on_submit(lambda text, w=which: self._on_color_bound(w, text))
        self.color_bounds = (low, high)
        self._color_hint = hint(0.79, y)
        y += row

        # A min and a max box per field, empty meaning "no bound".  Typing an
        # exact number is what this is actually used for -- a precision cut at
        # 20 nm, a z slab at +-150 -- and a slider cannot express that.  The
        # data range is shown beside each row so the numbers to type are known.
        # Laid out from the bottom up, so the fields read in their usual order.
        self.bounds: Dict[str, tuple] = {}
        self._hints: Dict[str, object] = {}
        for name in reversed(fields):
            low = TextBox(place(0.40, y, 0.17, box_h), name, initial="",
                          textalignment="right")
            high = TextBox(place(0.60, y, 0.17, box_h), "to", initial="",
                           textalignment="right")
            for box, which in ((low, 0), (high, 1)):
                box.on_submit(
                    lambda text, n=name, w=which: self._on_bound(n, w, text))
            self._hints[name] = hint(0.79, y)
            self.bounds[name] = (low, high)
            y += row
        for boxes in list(self.bounds.values()) + [self.color_bounds]:
            for box in boxes:
                box.label.set_fontsize(8)
                box.text_disp.set_fontsize(8)
        self._update_hints()

        if self.state.settings.color_field is not None:
            self._color_typed = list(
                self.state.settings.color_range
                or quantile_range(self.state.locs,
                                  self.state.settings.color_field, 0.01, 0.99))
            self._apply_color_range()

        self.image = None
        self._scalebar: List = []
        self._drag: Optional[Tuple[float, float]] = None
        self._render_seconds = 0.0      # how long the last render took
        self._rendered_at = 0.0
        self._timer = self.figure.canvas.new_timer(interval=self.REDRAW_DELAY_MS)
        self._timer.single_shot = True
        self._timer.add_callback(self._render_now)

        canvas = self.figure.canvas
        canvas.mpl_connect("scroll_event", self._on_scroll)
        canvas.mpl_connect("key_press_event", self._on_key)
        canvas.mpl_connect("button_press_event", self._on_press)
        canvas.mpl_connect("motion_notify_event", self._on_motion)
        canvas.mpl_connect("button_release_event", self._on_release)
        canvas.mpl_connect("resize_event", lambda _event: self._schedule())
        canvas.mpl_connect("close_event", lambda _event: self.close())

        self._apply_default_bounds()
        self._update_color_row()
        self.reset()

    # ------------------------------------------------------------------ set-up
    def _available_fields(self) -> List[str]:
        found = []
        for group in FILTER_FIELDS:
            for name in group:
                if name in self.state.locs:
                    found.append(name)
                    break
        return found

    def _available_color_fields(self) -> List[Tuple[str, Optional[str]]]:
        """Plain intensity, and every field this table could be coloured by."""
        found: List[Tuple[str, Optional[str]]] = [("intensity", None)]
        for label, names in COLOR_FIELDS:
            for name in names:
                if name in self.state.locs:
                    found.append((label, name))
                    break
        return found

    def _color_choice_for(self, field: Optional[str]) -> int:
        """Which radio entry a colour field is, adding one if it is not offered.

        `RenderSettings.color_field` takes any column, and the panel lists the
        few that are usually wanted; a table opened already coloured by
        something else gets that column as an entry of its own rather than a
        panel that disagrees with the image.
        """
        if field is None:
            return 0
        for i, (_, name) in enumerate(self.color_choices):
            if name == field:
                return i
        self.color_choices.append((field, field))
        return len(self.color_choices) - 1

    def _apply_default_bounds(self) -> None:
        """Put `DEFAULT_BOUNDS` on a table that has no bounds of its own yet.

        Once per table: a bound the user cleared must stay cleared, including
        across a switch to the grouped table and back.
        """
        current = self.state.current
        if getattr(current, "defaults_applied", False):
            return
        for field, (lo, hi) in DEFAULT_BOUNDS.items():
            if field in current.locs and field not in current.filter.ranges:
                current.filter.set(field, lo, hi)
        current.defaults_applied = True
        self._restore_bounds()

    def reset(self) -> None:
        """Show everything."""
        xrange, yrange = self.state.full_view()
        self._set_limits(xrange, (yrange[1], yrange[0]))
        self._render_now()

    # ----------------------------------------------------------------- drawing
    def _current_fov(self) -> FieldOfView:
        """A field of view matching the axes: one rendered pixel per screen pixel."""
        box = self.axes.get_window_extent()
        nx = max(int(round(box.width)), 16)
        ny = max(int(round(box.height)), 16)
        x0, x1 = self.axes.get_xlim()
        y1, y0 = self.axes.get_ylim()          # the y axis points down
        return FieldOfView.fit((x0, x1), (y0, y1), nx, ny)

    def _render_now(self) -> None:
        started = time.perf_counter()
        self.fov = self._current_fov()
        if len(self.state.locs) == 0:
            # A live view opens before the first block has been fitted.  The
            # limits are still set from the field of view, so the empty window
            # already frames what the first localizations will land in and the
            # image does not jump when they do.
            self._set_limits((self.fov.x0, self.fov.x1), (self.fov.y1, self.fov.y0))
            self._update_title()
            self.figure.canvas.draw_idle()
            return
        self._rendered = self.state.render(self.fov)
        rgb = self.state.display.apply(self._rendered)
        if self.image is None:
            self.image = self.axes.imshow(rgb, extent=self.fov.extent,
                                          origin="upper", interpolation="nearest")
            # imshow sets the aspect itself; put ours back
            self.axes.set_aspect("equal", adjustable="datalim")
        else:
            self.image.set_data(rgb)
            self.image.set_extent(self.fov.extent)
        self._set_limits((self.fov.x0, self.fov.x1), (self.fov.y1, self.fov.y0))
        self._draw_scalebar()
        self._update_title()
        self.figure.canvas.draw_idle()
        self._render_seconds = time.perf_counter() - started
        self._rendered_at = time.perf_counter()

    def _redisplay(self) -> None:
        """Re-apply the display settings to the planes we already have.

        The contrast is a quantile of the planes in hand, so it means the same
        thing at every zoom and needs no re-render to follow one.
        """
        if self._rendered is None:
            return
        self.image.set_data(self.state.display.apply(self._rendered))
        self.figure.canvas.draw_idle()

    def _draw_scalebar(self) -> None:
        for artist in self._scalebar:
            artist.remove()
        self._scalebar.clear()
        length = _nice_length(self.fov.x1 - self.fov.x0)
        pad = 0.03 * (self.fov.x1 - self.fov.x0)
        x = self.fov.x1 - pad - length
        y = self.fov.y1 - pad * (self.fov.y1 - self.fov.y0) / (self.fov.x1 - self.fov.x0)
        line, = self.axes.plot([x, x + length], [y, y], color="white", lw=3,
                               solid_capstyle="butt")
        unit = "nm" if "nm" in self.state.locs.metadata.get("units", "nm") else "px"
        label = self.axes.text(x + length / 2, y - 0.015 * (self.fov.y1 - self.fov.y0),
                               f"{length:g} {unit}", color="white", ha="center",
                               va="bottom", fontsize=9)
        self._scalebar += [line, label]

    def _update_title(self, note: str = "") -> None:
        shown = self._rendered.n_locs if self._rendered else 0
        kind = "grouped" if self.state.use_grouped else "localizations"
        if self.state.use_grouped and self.state.grouped_stale:
            kind += " (stale)"     # data has arrived since these were grouped
        title = note or (f"{shown:,} of {len(self.state.locs):,} {kind}   ·   "
                         f"{self.fov.pixelsize:.3g} per pixel   ·   "
                         f"{self.state.settings.mode}")
        if self.status:
            title += f"   ·   {self.status}"
        self.axes.set_title(title, fontsize=10)

    def _set_limits(self, xlim, ylim) -> None:
        self.axes.set_xlim(*xlim)
        self.axes.set_ylim(*ylim)
        self.axes.set_xticks([])
        self.axes.set_yticks([])

    def _schedule(self) -> None:
        """Re-render once the interaction has been still for a moment."""
        self._timer.stop()
        self._timer.start()

    def _refresh(self) -> None:
        """Follow a gesture: render inside it when that is fast, else defer.

        The throttle is the render's own cost, so a gesture never queues up
        renders it cannot keep up with -- at worst every other frame of the
        gesture is a rendered one.
        """
        if (self._render_seconds <= self.LIVE_RENDER_SECONDS
                and time.perf_counter() - self._rendered_at >= self._render_seconds):
            self._timer.stop()
            self._render_now()
        else:
            self.figure.canvas.draw_idle()
            self._schedule()

    # ------------------------------------------------------------------ events
    def _zoom(self, factor: float, cx: Optional[float] = None,
              cy: Optional[float] = None) -> None:
        x0, x1 = self.axes.get_xlim()
        y0, y1 = self.axes.get_ylim()
        cx = 0.5 * (x0 + x1) if cx is None else cx
        cy = 0.5 * (y0 + y1) if cy is None else cy
        self._set_limits((cx - (cx - x0) * factor, cx + (x1 - cx) * factor),
                         (cy - (cy - y0) * factor, cy + (y1 - cy) * factor))
        self._refresh()

    def _on_scroll(self, event) -> None:
        if event.inaxes is not self.axes:
            return
        self._zoom(1.2 ** -event.step, event.xdata, event.ydata)

    def _on_key(self, event) -> None:
        if event.key in ("+", "="):
            self._zoom(1 / 1.4)
        elif event.key in ("-", "_"):
            self._zoom(1.4)
        elif event.key == "r":
            self.reset()

    def _on_press(self, event) -> None:
        if event.inaxes is self.axes and event.button == 1:
            self._timer.stop()      # a queued render must not land mid-drag
            self._drag = (event.xdata, event.ydata)

    def _on_motion(self, event) -> None:
        if self._drag is None or event.xdata is None:
            return
        # the anchor is where the cursor grabbed the image, in data coordinates,
        # so panning is exact however far the view has moved
        dx, dy = self._drag[0] - event.xdata, self._drag[1] - event.ydata
        x0, x1 = self.axes.get_xlim()
        y0, y1 = self.axes.get_ylim()
        self._set_limits((x0 + dx, x1 + dx), (y0 + dy, y1 + dy))
        self._refresh()

    def _on_release(self, event) -> None:
        if self._drag is not None:
            self._drag = None
            self._schedule()

    def _on_bound(self, field: str, which: int, text: str) -> None:
        """One bound was typed.  Empty means no bound; junk is ignored."""
        try:
            value = float(text) if text.strip() else None
        except ValueError:
            self._restore_bounds()          # put the previous value back
            return
        lo, hi = self.state.filter.ranges.get(field, (None, None))
        lo, hi = (value, hi) if which == 0 else (lo, value)
        if lo is None and hi is None:
            self.state.filter.remove(field)
        else:
            try:
                self.state.filter.set(field, lo, hi)
            except KeyError:        # no such column (yet): a live view with no data
                self._restore_bounds()
                return
        self._schedule()

    def take_new_data(self, blocks: Sequence[Localizations]) -> int:
        """Append localizations and redraw.  Returns how many arrived.

        Nothing the user set is touched: the limits, the filter bounds, the
        contrast and the colour map are all left exactly as they are, so data
        appears inside the view being looked at instead of resetting it.  The
        one exception is the first block, which has to be framed because there
        was nothing to frame before.

        A redraw is deferred while the mouse is down, so an update can never
        interrupt a drag.
        """
        n = 0
        for block in blocks:
            n += self.state.append(block)
        if n == 0:
            return 0
        self._update_hints()
        self._update_color_row()    # the data range moves; what was typed does not
        if len(self.state.locs) == n:
            self._apply_default_bounds()    # the columns exist only now
        if self.state.use_grouped:
            self._update_title()    # the grouped table is stale, not different
            self.figure.canvas.draw_idle()
        elif len(self.state.locs) == n:
            self.reset()            # the first data: there was no view to keep
        elif self._drag is None:
            self._render_now()
        else:
            self._schedule()
        return n

    # ----------------------------------------------------------------- colour
    def _on_color(self, label: str) -> None:
        """Colour by a field, or by nothing.

        The LUT follows: an intensity image wants a brightness ramp and a coded
        field wants a full-hue one, and there is no reading a field off `hot`.
        The range starts at the field's own data range -- a z range means
        nothing for a frame number, so it is not carried across.
        """
        field = dict(self.color_choices).get(label)
        if field is not None and field not in self.state.locs:
            self._set_color_choice(self._color_label)   # not in this table
            self._update_color_row()
            return
        self._color_label = label
        self._set_color_choice(label)       # a no-op when it was clicked
        self.state.settings.color_field = field
        self.state.display.lut = INTENSITY_LUT if field is None else FIELD_LUT
        self._color_typed = ([None, None] if field is None
                             else list(quantile_range(self.state.locs, field,
                                                      0.01, 0.99)))
        self._apply_color_range()
        self._update_color_row()
        self._render_now()

    def _on_color_bound(self, which: int, text: str) -> None:
        """One end of the colour range was typed.  Empty means the data's own."""
        try:
            value = float(text) if text.strip() else None
        except ValueError:
            self._update_color_row()        # put the previous value back
            return
        self._color_typed[which] = value
        self._apply_color_range()
        self._schedule()

    def _apply_color_range(self) -> None:
        """Turn what was typed into the concrete range the renderer takes."""
        field = self.state.settings.color_field
        lo, hi = self._color_typed
        if field is None or (lo is None and hi is None):
            self.state.settings.color_range = None
            return
        values = np.asarray(self.state.locs[field], np.float64)
        self.state.settings.color_range = (
            float(np.nanmin(values)) if lo is None else lo,
            float(np.nanmax(values)) if hi is None else hi)

    def _update_color_row(self) -> None:
        """Show the field being coloured by, its range and the data's own."""
        field = self.state.settings.color_field
        for box, value in zip(self.color_bounds, self._color_typed):
            events, box.eventson = box.eventson, False
            try:
                box.set_val("" if value is None else f"{value:g}")
            finally:
                box.eventson = events
        if field is None:
            self._color_hint.set_text("(intensity)")
        elif field in self.state.locs and len(self.state.locs):
            lo, hi = quantile_range(self.state.locs, field, 0.01, 0.99)
            self._color_hint.set_text(f"{field}  [{lo:.4g} … {hi:.4g}]")
        else:
            self._color_hint.set_text(f"{field}  (not in this table)")

    def _set_color_choice(self, label: str) -> None:
        """Move the radio button without setting off its callback."""
        names = [name for name, _ in self.color_choices]
        if label not in names:
            return
        events, self.colors.eventson = self.colors.eventson, False
        try:
            self.colors.set_active(names.index(label))
        finally:
            self.colors.eventson = events

    def _on_contrast(self, value: float) -> None:
        self.state.display.contrast = float(value)
        self._redisplay()

    def _on_gamma(self, value: float) -> None:
        self.state.display.gamma = float(value)
        self._redisplay()

    def _on_switch(self, label: str) -> None:
        if label == "additive":
            self._on_additive(bool(self.switches.get_status()[1]))
        else:
            self._on_grouped()

    def _on_additive(self, on: bool) -> None:
        """Switch the composite for field-coloured images.

        ``sum`` is SMAP's: the colour planes are displayed directly, so a red
        and a cyan localization in one saturated pixel add to white.  ``hue``
        divides by the weight first, so that pixel stays a mid grey at full
        brightness and the hue never drifts towards white with density.  Below
        saturation the two are identical, so this only shows up where the image
        is brightest.  Both come from the same render.
        """
        self.state.display.color_mode = "sum" if on else "hue"
        self._redisplay()

    def _on_grouped(self) -> None:
        """Switch between the grouped and ungrouped tables.

        Each keeps its own filter, so the boxes are restored to what that table
        was last set to rather than carried across -- a precision cut that makes
        sense for single localizations is the wrong one for groups.
        """
        on = bool(self.switches.get_status()[0])
        if on and not self.state.has_grouped:      # the one slow path
            self._update_title("grouping...")
            self.figure.canvas.draw_idle()
            self.figure.canvas.flush_events()
        self.state.show_grouped(on, self.group_settings)
        self._apply_default_bounds()
        field = self.state.settings.color_field
        if field is not None and field not in self.state.locs:
            self._set_color_choice("intensity")
            self._on_color("intensity")     # renders; the rest is done there
            self._restore_bounds()
            return
        self._restore_bounds()
        self._render_now()

    def _update_hints(self) -> None:
        """Print the active table's data range beside each row.

        Grouping changes these -- a group's precision reaches well below any
        single localization's -- so they follow the switch.
        """
        for name, hint in self._hints.items():
            if name in self.state.locs:
                lo, hi = quantile_range(self.state.locs, name, 0.01, 0.99)
                hint.set_text(f"[{lo:.4g} … {hi:.4g}]")
            else:
                hint.set_text("(not in this table)")

    def _restore_bounds(self) -> None:
        """Show the active table's own filter in the boxes, silently."""
        self._update_hints()
        self._update_color_row()
        ranges = self.state.filter.ranges
        for field, boxes in self.bounds.items():
            for box, value in zip(boxes, ranges.get(field, (None, None))):
                events, box.eventson = box.eventson, False
                try:
                    box.set_val("" if value is None else f"{value:g}")
                finally:
                    box.eventson = events


def show(locs: Localizations, settings: Optional[RenderSettings] = None,
         display: Optional[DisplaySettings] = None, block: bool = True,
         group_settings: Optional[GroupSettings] = None) -> Viewer:
    """Open a viewer on a localization table."""
    import matplotlib.pyplot as plt

    viewer = Viewer(ViewState(locs, settings, display),
                    group_settings=group_settings)
    plt.show(block=block)
    return viewer
