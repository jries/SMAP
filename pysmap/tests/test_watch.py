"""Reading a file that is still being written.

The fake microscope here writes pages one at a time, sometimes into a second
file of the series, so the tests cover what actually goes wrong online: frames
that are not there yet, a page half-written, and an acquisition that just stops.
"""
import threading
import time

import numpy as np
import pytest
import tifffile

from smapfit.io.tiff import open_stack
from smapfit.io.watch import (WatchSettings, _read_new_pages,
                              open_growing_stack, watch_stack)


def _frame(i, shape=(16, 24)):
    """A frame that says which frame it is, so order can be checked."""
    return np.full(shape, i, dtype=np.uint16)


class Microscope:
    """Writes frames into an MM-style series, on demand or on a timer."""

    def __init__(self, directory, per_file=1000):
        self.directory = directory
        self.per_file = per_file
        self.written = 0

    def _path(self, n):
        stem = "acq_MMStack_Default"
        return self.directory / (f"{stem}.ome.tif" if n == 0
                                 else f"{stem}_{n}.ome.tif")

    def write(self, n=1):
        for _ in range(n):
            self._write_page(_frame(self.written))

    def _write_page(self, image):
        path = self._path(self.written // self.per_file)
        # ome=False, metadata=None: tifffile refuses to append to a file
        # carrying its own metadata, and writing page by page is the point
        with tifffile.TiffWriter(path, append=path.exists(), ome=False) as tif:
            tif.write(image, metadata=None)
        self.written += 1

    def run(self, n, interval, delay=0.0):
        """Write ``n`` frames in the background, one every ``interval``."""
        def writer():
            time.sleep(delay)
            for _ in range(n):
                self.write()
                time.sleep(interval)
        thread = threading.Thread(target=writer, daemon=True)
        thread.start()
        return thread


def _collect(stream):
    """Drain a frame stream into (first_frame, frame_number) pairs."""
    seen = []
    for first, block in stream:
        assert block.ndim == 3
        for i, image in enumerate(block):
            seen.append((first + i, int(image.flat[0])))
    return seen


FAST = WatchSettings(poll=0.02, timeout=0.4, appear_timeout=2.0)


def test_frames_written_before_the_watch_starts_are_read(tmp_path):
    scope = Microscope(tmp_path)
    scope.write(7)
    seen = _collect(watch_stack(tmp_path, chunk=3, settings=FAST))
    assert seen == [(i, i) for i in range(7)]


def test_frames_arriving_during_the_watch_are_read_in_order(tmp_path):
    scope = Microscope(tmp_path)
    scope.write(2)
    scope.run(20, interval=0.01, delay=0.05)
    seen = _collect(watch_stack(tmp_path, chunk=4, settings=FAST))
    assert [n for _, n in seen] == list(range(22))
    assert [i for i, _ in seen] == list(range(22))   # frame indices are continuous


def test_the_watch_follows_the_series_into_the_next_file(tmp_path):
    scope = Microscope(tmp_path, per_file=5)
    scope.write(12)
    assert (tmp_path / "acq_MMStack_Default_2.ome.tif").exists()
    seen = _collect(watch_stack(tmp_path, chunk=5, settings=FAST))
    assert [n for _, n in seen] == list(range(12))


def test_the_watch_ends_when_the_acquisition_stops(tmp_path):
    scope = Microscope(tmp_path)
    scope.write(3)
    t0 = time.monotonic()
    seen = _collect(watch_stack(tmp_path, chunk=100, settings=FAST))
    elapsed = time.monotonic() - t0
    assert [n for _, n in seen] == [0, 1, 2]
    assert FAST.timeout <= elapsed < FAST.timeout + 2.0


def test_a_partial_block_is_handed_over_rather_than_held_back(tmp_path):
    """A sparse acquisition must not wait for a full block to appear."""
    scope = Microscope(tmp_path)
    scope.write(3)
    scope.run(2, interval=0.15, delay=0.05)
    blocks = [len(block) for _, block in watch_stack(tmp_path, chunk=1000,
                                                     settings=FAST)]
    assert sum(blocks) == 5
    assert len(blocks) > 1


def test_the_page_being_written_is_held_back(tmp_path):
    """The newest page of a growing file may be half there, so it waits.

    Held back, not dropped: once the acquisition has stopped, nothing can
    follow that page and the watch reads it on its way out -- which is what the
    other tests, all of which get every frame, check from the outside.
    """
    scope = Microscope(tmp_path)
    scope.write(4)
    path = tmp_path / "acq_MMStack_Default.ome.tif"

    images, page_i, exhausted = _read_new_pages(path, 0, hold_last=True)
    assert [int(im.flat[0]) for im in images] == [0, 1, 2]
    assert page_i == 3 and exhausted

    images, page_i, _ = _read_new_pages(path, page_i, hold_last=False)
    assert [int(im.flat[0]) for im in images] == [3]
    assert page_i == 4


def test_a_file_caught_mid_write_is_a_not_yet_rather_than_a_failure(tmp_path):
    broken = tmp_path / "half.ome.tif"
    broken.write_bytes(b"II*\x00" + b"\x00" * 40)      # a header and nothing else
    images, page_i, _ = _read_new_pages(broken, 0, hold_last=False)
    assert (images, page_i) == ([], 0)          # no exception escapes; try again

    images, page_i, _ = _read_new_pages(tmp_path / "gone.ome.tif", 0,
                                        hold_last=False)
    assert (images, page_i) == ([], 0)


def test_waiting_for_a_file_that_does_not_exist_times_out(tmp_path):
    settings = WatchSettings(poll=0.02, timeout=0.1, appear_timeout=0.2)
    assert _collect(watch_stack(tmp_path / "nothing", settings=settings)) == []
    with pytest.raises(TimeoutError):
        open_growing_stack(tmp_path / "nothing", settings=settings)


def test_a_file_that_appears_late_is_picked_up(tmp_path):
    scope = Microscope(tmp_path)
    scope.run(4, interval=0.02, delay=0.2)
    seen = _collect(watch_stack(tmp_path, chunk=2, settings=FAST))
    assert [n for _, n in seen] == [0, 1, 2, 3]


def test_stopping_early_ends_the_stream(tmp_path):
    scope = Microscope(tmp_path)
    scope.write(4)
    stop = threading.Event()
    stream = watch_stack(tmp_path, chunk=1, settings=FAST, stop_event=stop)
    assert next(stream)[0] == 0
    stop.set()
    assert list(stream) == []


def test_start_and_stop_select_a_range(tmp_path):
    scope = Microscope(tmp_path)
    scope.write(10)
    seen = _collect(watch_stack(tmp_path, chunk=3, settings=FAST,
                                start=2, stop=7))
    assert [n for _, n in seen] == [2, 3, 4, 5, 6]


def test_a_growing_stack_can_be_opened_for_its_metadata(tmp_path):
    scope = Microscope(tmp_path)
    scope.run(3, interval=0.02, delay=0.1)
    source = open_growing_stack(tmp_path, settings=FAST)
    assert source.shape == (16, 24)
    assert source.dtype == np.uint16
    assert source.n_frames >= 1        # only what was there when it opened
