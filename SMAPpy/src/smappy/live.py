"""Fitting an acquisition while it runs, with the viewer open on the result.

The two halves already existed and are deliberately not merged here: the
fitter is `LocalizationEngine` driven by a file that is still being written
(`smappy.io.watch`), the view is `ViewState`, which can be appended to.  All
this module adds is the join between them, which is a **queue and a timer**.

* The fitter runs in its own thread and puts each finished block on a queue.
  It never touches the view.  Nothing needs a lock, because a block of
  localizations is handed over and never looked at again by the thread that
  made it -- and because the C++ that dominates the fit releases the GIL, that
  thread does not compete with the window for it either.
* The viewer drains the queue on a timer, every few seconds, and appends.
  Between two of those moments the window belongs entirely to the user.

Nothing an update does changes a setting.  Filter bounds, zoom, contrast, the
colour map and the choice of table all survive an append untouched, because
appending extends the table, the filter masks and the index rather than
rebuilding any of them (see `ViewState.append`).

An HDF5 file is written at the same time and is the actual result; the view is
a look at what is going into it, not the record.  If the viewer is closed the
fit keeps writing until the acquisition ends, unless it is stopped too.
"""

from __future__ import annotations

import queue
import sys
import threading
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .io.watch import WatchSettings, open_growing_stack
from .locs import Localizations
from .metadata import CameraMetadata
from .pipeline import FitSettings, LocalizationEngine, prefetch, provenance
from .render import DisplaySettings, RenderSettings
from .viewer import COLOR_FIELDS, FILTER_FIELDS, ViewState, Viewer


@dataclass
class LiveSettings:
    """The timing of an online run: how often things are handed on."""

    chunk: int = 100             # frames read from the file at a time
    update_seconds: float = 3.0  # how often the view takes in what is queued
    flush_seconds: float = 5.0   # fit the buffered ROIs at least this often
    read_ahead: int = 2          # blocks the reader may run ahead of the fit
    watch: WatchSettings = field(default_factory=WatchSettings)


class LiveFit:
    """The pipeline, running in a thread, handing finished blocks over a queue.

    The loop is `fit_stack`'s with one addition: the engine is flushed on a
    timer as well as when its ROI buffer fills.  Without that, a sparse
    acquisition would show nothing for minutes -- the buffer holds 15000 ROIs
    and a quiet sample does not produce them quickly.
    """

    def __init__(self, frames: Iterable[Tuple[int, np.ndarray]],
                 camera: CameraMetadata, finder, model,
                 settings: Optional[FitSettings] = None,
                 live: Optional[LiveSettings] = None,
                 sink=None, stop_event: Optional[threading.Event] = None):
        self.engine = LocalizationEngine(camera, finder, model,
                                         settings or FitSettings())
        self.frames = frames
        self.live = live or LiveSettings()
        self.sink = sink
        self.stop_event = stop_event or threading.Event()
        self.results: "queue.Queue[Localizations]" = queue.Queue()
        self.finished = threading.Event()
        self.error: Optional[BaseException] = None
        self.n_emitted = 0
        self._thread: Optional[threading.Thread] = None
        self._t0 = time.monotonic()

    # ----------------------------------------------------------------- running
    def start(self) -> "LiveFit":
        self._t0 = time.monotonic()
        self._thread = threading.Thread(target=self._run, name="smappy-live",
                                        daemon=True)
        self._thread.start()
        return self

    def _run(self) -> None:
        try:
            last_flush = time.monotonic()
            for first_frame, block in prefetch(self.frames,
                                               self.live.read_ahead):
                if self.stop_event.is_set():
                    break
                self._emit(self.engine.push(block, first_frame))
                if time.monotonic() - last_flush >= self.live.flush_seconds:
                    self._emit(self.engine.flush())
                    last_flush = time.monotonic()
            self._emit(self.engine.flush())
        except BaseException as error:       # reported, not swallowed
            self.error = error
            traceback.print_exc(file=sys.stderr)
        finally:
            self.finished.set()

    def _emit(self, locs: Optional[Localizations]) -> None:
        if locs is None or len(locs) == 0:
            return
        if self.sink is not None:
            self.sink(locs)                  # the file first: it is the result
        self.n_emitted += len(locs)
        self.results.put(locs)

    def stop(self, timeout: float = 5.0) -> None:
        """Ask the fit to end and wait for it, at most ``timeout``."""
        self.stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout)

    # ----------------------------------------------------------------- reading
    def drain(self) -> List[Localizations]:
        """Every block queued since the last call."""
        blocks = []
        while True:
            try:
                blocks.append(self.results.get_nowait())
            except queue.Empty:
                return blocks

    @property
    def running(self) -> bool:
        return not self.finished.is_set()

    @property
    def status(self) -> str:
        """One line for the title: how far the fit has got."""
        stats = self.engine.stats
        elapsed = max(time.monotonic() - self._t0, 1e-9)
        state = "fitting" if self.running else "acquisition ended"
        if self.error is not None:
            state = f"stopped: {type(self.error).__name__}"
        return (f"{state} · {stats['frames']:,} frames "
                f"({stats['frames'] / elapsed:.0f}/s)")


class LiveViewer(Viewer):
    """A viewer that takes in what a `LiveFit` has queued, on a timer."""

    def __init__(self, state: ViewState, fit: LiveFit,
                 live: Optional[LiveSettings] = None, on_close=None, **kwargs):
        super().__init__(state, **kwargs)
        self.fit = fit
        self.live = live or LiveSettings()
        self._on_close = on_close
        self._closed = False
        self._updates = self.figure.canvas.new_timer(
            interval=max(int(self.live.update_seconds * 1000), 100))
        self._updates.add_callback(self.update)
        self._updates.start()
        self.figure.canvas.mpl_connect("close_event", lambda _e: self.close())

    def update(self) -> int:
        """One tick: take what is queued, redraw, and keep the status current."""
        was_running = self.fit.running
        blocks = self.fit.drain()
        self.status = self.fit.status
        n = self.take_new_data(blocks)
        if n == 0 and self.fov is not None:
            self._update_title()             # the status still moves
            self.figure.canvas.draw_idle()
        if not was_running:                  # drained after it finished: done
            self._updates.stop()
        return n

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._updates.stop()
        if self._on_close is not None:
            self._on_close()


# ------------------------------------------------------------------ putting it up
def live_view(data, camera: CameraMetadata, finder, model,
              settings: Optional[FitSettings] = None,
              output=None,
              live: Optional[LiveSettings] = None,
              render_settings: Optional[RenderSettings] = None,
              display: Optional[DisplaySettings] = None,
              group_settings=None, extent: Optional[Sequence[float]] = None,
              stop_fit_on_close: bool = True, block: bool = True) -> LiveViewer:
    """Fit an acquisition as it is written and watch the image build up.

    ``data`` is the growing Micro-Manager TIFF (or the directory holding it),
    ``output`` the HDF5 to write -- the result; the view is a look at it.  The
    call waits for the first frames to appear, opens the window, and returns
    when it is closed (``block=False`` returns at once, for a notebook).

    The field of view the image is framed in comes from the camera, so it is
    the same from the first localization on and the picture does not rescale
    under the user as data arrives; pass ``extent`` to override it.
    """
    import matplotlib.pyplot as plt

    from .io.hdf5 import LocalizationWriter

    live = live or LiveSettings()
    settings = settings or FitSettings()
    stop = threading.Event()

    source = data
    if isinstance(data, (str, Path)):
        source = open_growing_stack(data, live.watch, stop_event=stop)
    frames = source.watch(chunk=live.chunk, settings=live.watch, stop_event=stop)

    writer = None
    if output is not None:
        writer = LocalizationWriter(output)
        writer.set_metadata(provenance(camera, finder, model, settings,
                                       source=getattr(source, "files", [data])[0]))

    fit = LiveFit(frames, camera, finder, model, settings, live,
                  sink=None if writer is None else writer.append,
                  stop_event=stop)

    if extent is None:
        extent = camera_extent(camera, source.shape, settings.output_unit)
    state = ViewState(Localizations(), render_settings, display,
                      live=True, extent=extent)

    def on_close():
        if stop_fit_on_close:
            fit.stop()
        if writer is not None and not fit.running:
            writer.close()

    viewer = LiveViewer(state, fit, live, on_close=on_close,
                        filter_fields=expected_filter_fields(settings.output_unit),
                        color_fields=expected_color_fields(settings.output_unit),
                        group_settings=group_settings)
    fit.start()
    plt.show(block=block)
    if block:
        viewer.close()                # stops the fit, if that is what was asked
        if not stop_fit_on_close:
            fit.finished.wait()       # closing the window only closes the window
        fit.stop()
        if writer is not None:
            writer.close()
    return viewer


def expected_filter_fields(output_unit: str) -> List[str]:
    """The filter rows to put up before the first block says what the columns are.

    A live window opens on an empty table, so the fields cannot be read off it;
    they follow from the unit the fit will write.  A row for a column this fit
    does not produce -- z without a 3D calibration -- says so beside itself
    rather than being guessed away.
    """
    suffix = "_pix" if output_unit == "pixel" else "_nm"
    fields = []
    for group in FILTER_FIELDS:
        match = [name for name in group if name.endswith(suffix)]
        fields.append(match[0] if match else group[0])
    return fields


def expected_color_fields(output_unit: str) -> List[Tuple[str, Optional[str]]]:
    """The "colour by" choices to offer before the first block has arrived.

    Same reasoning as `expected_filter_fields`: a live window opens on an empty
    table, so the columns have to follow from the unit the fit will write.
    Choosing one the fit does not produce says so and moves back.
    """
    suffix = "_pix" if output_unit == "pixel" else "_nm"
    choices: List[Tuple[str, Optional[str]]] = [("intensity", None)]
    for label, names in COLOR_FIELDS:
        match = [name for name in names if name.endswith(suffix)]
        choices.append((label, match[0] if match else names[0]))
    return choices


def camera_extent(camera: CameraMetadata, shape, output_unit: str = "nm"
                  ) -> Optional[Tuple[float, float, float, float]]:
    """The rectangle the camera can see, in the units the table will use.

    Fixing the view to this instead of to the data means the image is framed
    the same way from the first localization to the last.
    """
    height, width = int(shape[0]), int(shape[1])
    x0, y0 = camera.roi_offset
    if camera.roi is not None and len(camera.roi) >= 4:
        width, height = int(camera.roi[2]), int(camera.roi[3])
    scale = 1.0
    if output_unit != "pixel":
        if not camera.pixelsize_um:
            return None
        scale = float(camera.pixelsize_um) * 1000.0
    return (x0 * scale, y0 * scale, (x0 + width) * scale, (y0 + height) * scale)
