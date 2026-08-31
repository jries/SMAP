"""The join between a running fit and an open viewer.

The end-to-end test is the point: a fake microscope writes a stack frame by
frame, the fitter follows it, and the view is checked to hold everything -- and
to hold nothing the user set, which is the property that is easy to break.
"""
import threading
import time

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from smappy.detect import AbsoluteCutoff, DoGFilter, PeakFinder  # noqa: E402
from smappy.io.hdf5 import load_localizations                    # noqa: E402
from smappy.io.watch import WatchSettings                        # noqa: E402
from smappy.live import (LiveFit, LiveSettings, LiveViewer,      # noqa: E402
                          camera_extent, expected_filter_fields)
from smappy.locs import Localizations                            # noqa: E402
from smappy.metadata import CameraMetadata                       # noqa: E402
from smappy.pipeline import FitSettings                          # noqa: E402
from smappy.psf import GaussianPSF                               # noqa: E402
from smappy.viewer import ViewState                              # noqa: E402

from test_watch import Microscope                                 # noqa: E402

SHAPE = (48, 48)
CAMERA = CameraMetadata(conversion=1.0, offset=100.0, pixelsize_um=0.1,
                        em_on=False, emgain=1.0)
FAST = LiveSettings(chunk=4, update_seconds=0.05, flush_seconds=0.05,
                    watch=WatchSettings(poll=0.02, timeout=0.3,
                                        appear_timeout=2.0))


class Blinking(Microscope):
    """Frames with one bright emitter each, so every frame yields a fit."""

    def write(self, n=1):
        rng = np.random.default_rng(self.written)
        for _ in range(n):
            y, x = np.mgrid[0:SHAPE[0], 0:SHAPE[1]]
            cx, cy = rng.uniform(8, SHAPE[1] - 8), rng.uniform(8, SHAPE[0] - 8)
            frame = 100.0 + 2000.0 * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / 2.0)
            image = rng.poisson(np.maximum(frame - 100.0, 0)).astype(np.float64) + 100.0
            self._write_page(image.astype(np.uint16))


def _pipeline():
    finder = PeakFinder(DoGFilter(1.2), AbsoluteCutoff(40.0))
    return finder, GaussianPSF(sigma=1.2)


def _run(directory, out=None, frames=12, live=FAST):
    from smappy.io.watch import watch_stack
    finder, model = _pipeline()
    stream = watch_stack(directory, chunk=live.chunk, settings=live.watch)
    fit = LiveFit(stream, CAMERA, finder, model,
                  FitSettings(roisize=9, output_unit="nm"), live)
    return fit


def test_the_fit_follows_a_stack_that_is_still_being_written(tmp_path):
    scope = Blinking(tmp_path)
    scope.write(4)
    scope.run(12, interval=0.02, delay=0.05)

    fit = _run(tmp_path).start()
    collected = Localizations()
    while fit.running or not fit.results.empty():
        for block in fit.drain():
            collected.extend(block)
        time.sleep(0.02)
    for block in fit.drain():
        collected.extend(block)

    assert fit.error is None
    assert fit.engine.stats["frames"] == 16
    assert len(collected) > 8                       # most frames give a fit
    assert collected["frame"].max() < 16
    assert np.all(np.diff(collected["frame"]) >= 0)  # blocks arrive in order


def test_blocks_are_flushed_on_a_timer_so_sparse_data_still_shows(tmp_path):
    """The ROI buffer holds thousands; a quiet sample must not wait for it."""
    scope = Blinking(tmp_path)
    scope.write(3)
    fit = _run(tmp_path).start()
    fit.finished.wait(timeout=10)
    assert fit.n_emitted > 0
    assert fit.engine.stats["rois"] < fit.engine.settings.block_rois()


def test_the_viewer_takes_in_new_data_without_changing_a_setting(tmp_path):
    scope = Blinking(tmp_path)
    scope.write(6)
    fit = _run(tmp_path).start()
    fit.finished.wait(timeout=10)
    assert fit.error is None

    extent = camera_extent(CAMERA, SHAPE, "nm")
    state = ViewState(Localizations(), live=True, extent=extent)
    viewer = LiveViewer(state, fit, FAST,
                        filter_fields=expected_filter_fields("nm"))
    try:
        assert viewer.update() > 0                  # the first block frames it
        view = viewer.axes.get_xlim(), viewer.axes.get_ylim()
        viewer.state.filter.set("loc_precision_nm", None, 40.0)
        viewer.state.display.gamma = 0.5
        before = len(state.locs)

        state.append(_more_locs(state.locs, 500))
        viewer.take_new_data([])                    # nothing queued: a no-op tick

        assert len(state.locs) == before + 500
        assert (viewer.axes.get_xlim(), viewer.axes.get_ylim()) == view
        assert state.filter.ranges["loc_precision_nm"] == (None, 40.0)
        assert state.display.gamma == 0.5
        assert "fitting" in viewer.status or "ended" in viewer.status
    finally:
        viewer.close()


def _more_locs(like: Localizations, n: int) -> Localizations:
    rng = np.random.default_rng(7)
    return Localizations({name: rng.choice(np.asarray(like[name]), n)
                          for name in like.keys()}, dict(like.metadata))


def test_the_view_is_framed_by_the_camera_before_any_data_arrives(tmp_path):
    extent = camera_extent(CAMERA, SHAPE, "nm")
    state = ViewState(Localizations(), live=True, extent=extent)
    empty = state.full_view()
    assert empty[0][0] < 0.0 and empty[0][1] > SHAPE[1] * 100.0 * 0.9

    scope = Blinking(tmp_path)
    scope.write(6)
    fit = _run(tmp_path).start()
    fit.finished.wait(timeout=10)
    for block in fit.drain():
        state.append(block)
    assert len(state.locs) > 0
    assert state.full_view() == empty          # the frame did not move


def test_the_result_file_and_the_view_agree(tmp_path):
    from smappy.io.hdf5 import LocalizationWriter

    scope = Blinking(tmp_path)
    scope.write(8)
    out = tmp_path / "out.h5"
    with LocalizationWriter(out) as writer:
        fit = _run(tmp_path)
        fit.sink = writer.append
        fit.start()
        fit.finished.wait(timeout=10)
        state = ViewState(Localizations(), live=True)
        for block in fit.drain():
            state.append(block)

    written = load_localizations(out)
    assert len(written) == len(state.locs) == fit.n_emitted
    assert np.allclose(written["x_nm"], state.locs["x_nm"])


def test_stopping_ends_the_fit(tmp_path):
    scope = Blinking(tmp_path)
    scope.write(2)
    scope.run(200, interval=0.05, delay=0.0)
    fit = _run(tmp_path, live=LiveSettings(
        chunk=1, watch=WatchSettings(poll=0.02, timeout=30.0,
                                     appear_timeout=2.0))).start()
    time.sleep(0.3)
    t0 = time.monotonic()
    fit.stop(timeout=5.0)
    assert time.monotonic() - t0 < 5.0
    assert not fit.running


def test_live_view_wires_the_whole_thing_up(tmp_path):
    """The entry point itself: an empty window, then a fit filling it in."""
    import matplotlib.pyplot as plt

    from smappy.live import live_view

    data = tmp_path / "acq"
    data.mkdir()
    scope = Blinking(data)
    scope.write(2)
    scope.run(10, interval=0.02, delay=0.05)

    out = tmp_path / "live.h5"
    finder, model = _pipeline()
    viewer = live_view(data, CAMERA, finder, model,
                       FitSettings(roisize=9, output_unit="nm"),
                       output=out, live=FAST, block=False)
    try:
        assert len(viewer.state.locs) == 0          # the window opens empty
        assert viewer.axes.get_xlim()[1] > SHAPE[1] * 100.0 * 0.9   # framed already
        # the first draw settles the axes box (equal aspect in a new window),
        # so re-frame once against it before recording what must not move
        viewer.figure.canvas.draw()
        viewer._render_now()
        framed = viewer.axes.get_xlim(), viewer.axes.get_ylim()

        deadline = time.monotonic() + 20
        while viewer.fit.running and time.monotonic() < deadline:
            viewer.update()
            time.sleep(0.05)
        viewer.update()

        assert viewer.fit.error is None
        assert len(viewer.state.locs) > 0
        assert (viewer.axes.get_xlim(), viewer.axes.get_ylim()) == framed
        assert "loc_precision_nm" in viewer.bounds   # the rows are up from the start
    finally:
        viewer.close()
        viewer.fit.stop()
        plt.close(viewer.figure)

    written = load_localizations(out)
    assert len(written) == len(viewer.state.locs) == viewer.fit.n_emitted
    assert written.metadata["fit"]["roisize"] == 9
