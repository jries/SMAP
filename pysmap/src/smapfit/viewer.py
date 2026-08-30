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
slab), and the useful ranges span decades.  An empty box means no bound, so
opening a file never silently drops localizations.  The data range is printed
beside each row.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .filter import LocFilter, quantile_range
from .group import GroupSettings, group
from .locs import Localizations
from .render import (DisplaySettings, FieldOfView, RenderSettings, RenderedImage,
                     positions, render_locs)
from .spatial import SpatialIndex

# The fields worth a filter control, best alternative first.  Everything else
# stays available through `LocFilter` without cluttering the window.
FILTER_FIELDS: Sequence[Tuple[str, ...]] = (
    ("loc_precision_nm", "loc_precision_pix"),
    ("z_nm",),
    ("sigma_nm", "sigma_pix"),
    ("logl_rel",),
    ("frame",),
)


class LocSet:
    """One table with the things that are derived from it and nothing else.

    Grouped and ungrouped localizations are different tables, so each needs its
    own filter and its own spatial index.  Keeping both alive is the point:
    linking costs a third of a second and the filter caches would have to be
    rebuilt on every switch otherwise.
    """

    def __init__(self, locs: Localizations, filter: Optional[LocFilter] = None):
        self.locs = locs
        self.filter = filter if filter is not None else LocFilter(locs)
        x, y = positions(locs)
        self.index = SpatialIndex(x, y)

        # for the query margin: how far the widest blob reaches beyond the view
        self.median_precision = 0.0
        for name in ("loc_precision_nm", "loc_precision_pix"):
            if name in locs:
                values = np.asarray(locs[name], np.float64)
                values = values[np.isfinite(values)]
                if values.size:
                    self.median_precision = float(np.median(values))
                break

    def __len__(self) -> int:
        return len(self.locs)


class ViewState:
    """Everything needed to render a view, without a window around it."""

    def __init__(self, locs: Localizations,
                 settings: Optional[RenderSettings] = None,
                 display: Optional[DisplaySettings] = None,
                 filter: Optional[LocFilter] = None,
                 grouped: Optional[Localizations] = None):
        self.settings = settings or RenderSettings()
        self.display = display or DisplaySettings()
        self.n_threads = 0
        self.sets: Dict[str, LocSet] = {"ungrouped": LocSet(locs, filter)}
        if grouped is not None:
            self.sets["grouped"] = LocSet(grouped)
        self.use_grouped = False

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
        """Build the grouped table, once.  This is the slow part."""
        if "grouped" not in self.sets:
            grouped, _ = group(self.sets["ungrouped"].locs, settings)
            self.sets["grouped"] = LocSet(grouped)
        return self.sets["grouped"]

    def show_grouped(self, on: bool,
                     settings: Optional[GroupSettings] = None) -> None:
        if on:
            self.group(settings)
        self.use_grouped = bool(on)

    def full_view(self, margin_fraction: float = 0.01) -> Tuple[tuple, tuple]:
        """The coordinate ranges covering every localization."""
        ix = self.index
        x1 = ix.x0 + ix.n_cols * ix.cell_size
        y1 = ix.y0 + ix.n_rows * ix.cell_size
        mx, my = (x1 - ix.x0) * margin_fraction, (y1 - ix.y0) * margin_fraction
        return (ix.x0 - mx, x1 + mx), (ix.y0 - my, y1 + my)

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


class Viewer:
    """A matplotlib window over a :class:`ViewState`.

    Scroll or pinch to zoom about the cursor, ``+``/``-`` to zoom about the
    centre, drag to pan, ``r`` to reset.  Each filter field takes a minimum and
    a maximum, applied on Enter; an empty box means no bound.  The two sliders
    set the contrast and gamma, which re-display without re-rendering.
    "grouped" switches to the grouped table, building it on first use;
    "additive" switches the field-colour composite to SMAP's.
    """

    REDRAW_DELAY_MS = 150

    def __init__(self, state: ViewState, figsize=(9.5, 8.5),
                 filter_fields: Optional[Sequence[str]] = None,
                 group_settings: Optional[GroupSettings] = None):
        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, Slider, TextBox

        self.state = state
        self.group_settings = group_settings
        self.fov: Optional[FieldOfView] = None
        self._rendered: Optional[RenderedImage] = None

        fields = list(filter_fields) if filter_fields is not None \
            else self._available_fields()

        self.figure = plt.figure(figsize=figsize)
        rows = len(fields) + 2                      # + contrast and gamma
        bottom = 0.06 + 0.045 * rows
        self.axes = self.figure.add_axes([0.08, bottom, 0.90, 0.94 - bottom])
        self.axes.set_facecolor("black")
        self.axes.set_aspect("equal")

        # A min and a max box per field, empty meaning "no bound".  Typing an
        # exact number is what this is actually used for -- a precision cut at
        # 20 nm, a z slab at +-150 -- and a slider cannot express that.  The
        # data range is shown beside each row so the numbers to type are known.
        self.bounds: Dict[str, tuple] = {}
        self._hints: Dict[str, object] = {}
        for i, name in enumerate(fields):
            y = 0.055 + 0.045 * (rows - 1 - i)
            low = TextBox(self.figure.add_axes([0.30, y, 0.09, 0.032]), name,
                          initial="", textalignment="right")
            high = TextBox(self.figure.add_axes([0.42, y, 0.09, 0.032]), "to",
                           initial="", textalignment="right")
            for box, which in ((low, 0), (high, 1)):
                box.on_submit(lambda text, n=name, w=which: self._on_bound(n, w, text))
            self._hints[name] = self.figure.text(0.53, y + 0.008, "", fontsize=8,
                                                 color="0.45")
            self.bounds[name] = (low, high)
        self._update_hints()

        # SMAP's dynamic contrast: saturate 10^-p of the pixels.  The useful
        # range spans decades, which is why the slider moves the exponent.
        ax = self.figure.add_axes([0.20, 0.055 + 0.045, 0.62, 0.03])
        self.contrast = Slider(ax, "contrast (p)", 1.0, 5.0,
                               valinit=state.display.contrast)
        self.contrast.on_changed(self._on_contrast)
        ax = self.figure.add_axes([0.20, 0.055, 0.62, 0.03])
        self.gamma = Slider(ax, "gamma", 0.2, 2.0, valinit=state.display.gamma)
        self.gamma.on_changed(self._on_gamma)

        # "grouped" rebuilds the table (once); "additive" only changes how the
        # colour planes are composited, so it re-displays without re-rendering
        ax = self.figure.add_axes([0.845, bottom - 0.085, 0.145, 0.075])
        ax.set_frame_on(False)
        self.switches = CheckButtons(
            ax, ["grouped", "additive"],
            [state.use_grouped, state.display.color_mode == "sum"])
        self.switches.on_clicked(self._on_switch)
        self.grouped = self.switches   # the switch panel, kept under both names

        self.image = None
        self._scalebar: List = []
        self._drag: Optional[Tuple[float, float]] = None
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
        self.fov = self._current_fov()
        self._rendered = self.state.render(self.fov)
        rgb = self.state.display.apply(self._rendered)
        if self.image is None:
            self.image = self.axes.imshow(rgb, extent=self.fov.extent,
                                          origin="upper", interpolation="nearest")
        else:
            self.image.set_data(rgb)
            self.image.set_extent(self.fov.extent)
        self._set_limits((self.fov.x0, self.fov.x1), (self.fov.y1, self.fov.y0))
        self._draw_scalebar()
        self._update_title()
        self.figure.canvas.draw_idle()

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
        self.axes.set_title(
            note or f"{shown:,} of {len(self.state.locs):,} {kind}   ·   "
                    f"{self.fov.pixelsize:.3g} per pixel   ·   "
                    f"{self.state.settings.mode}", fontsize=10)

    def _set_limits(self, xlim, ylim) -> None:
        self.axes.set_xlim(*xlim)
        self.axes.set_ylim(*ylim)
        self.axes.set_xticks([])
        self.axes.set_yticks([])

    def _schedule(self) -> None:
        """Re-render once the interaction has been still for a moment."""
        self._timer.stop()
        self._timer.start()

    # ------------------------------------------------------------------ events
    def _zoom(self, factor: float, cx: Optional[float] = None,
              cy: Optional[float] = None) -> None:
        x0, x1 = self.axes.get_xlim()
        y0, y1 = self.axes.get_ylim()
        cx = 0.5 * (x0 + x1) if cx is None else cx
        cy = 0.5 * (y0 + y1) if cy is None else cy
        self._set_limits((cx - (cx - x0) * factor, cx + (x1 - cx) * factor),
                         (cy - (cy - y0) * factor, cy + (y1 - cy) * factor))
        self.figure.canvas.draw_idle()
        self._schedule()

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
        self.figure.canvas.draw_idle()

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
            self.state.filter.set(field, lo, hi)
        self._schedule()

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
