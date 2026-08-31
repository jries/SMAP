"""Reading an acquisition while the microscope is still writing it.

`ImageSource.frames` walks the pages a file had when it was opened, which is
everything an offline analysis needs and nothing an online one does.  Here the
file is re-examined instead: read what has arrived, hand it on, wait, look
again, and stop once nothing new has appeared for a while -- which is how an
acquisition ends, since Micro-Manager does not announce it and the frame count
in the summary metadata is what was *planned*, not what was written.

Three things make this safe against a writer working on the same file:

* **The last page of the file being written is held back** until another page
  follows it, so a page whose pixels are still being written is never read.
* **Every read is retried.**  A file caught mid-write raises rather than
  returning half an image; that is a "not yet", not a failure.
* **Files are re-listed on every poll**, so the ``_1.ome.tif`` that appears
  when the current one fills up is picked up without reopening anything.

The result is the same ``(first_frame, block)`` stream `ImageSource.frames`
yields, so the pipeline does not know the difference.
"""

from __future__ import annotations

import contextlib
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import numpy as np
import tifffile

from .tiff import ImageSource, _first_of_series, _series_files, open_stack


@dataclass
class WatchSettings:
    """When to look again, and when to give up."""

    poll: float = 1.0        # seconds between looks at the file
    timeout: float = 30.0    # stop after this long with no new frame
    appear_timeout: float = 60.0   # how long to wait for the file to show up
    hold_last: bool = True   # never read the newest page of a growing file


def watch_stack(path, chunk: int = 100, settings: Optional[WatchSettings] = None,
                start: int = 0, stop: Optional[int] = None,
                stop_event=None, on_wait=None
                ) -> Iterator[Tuple[int, np.ndarray]]:
    """Yield ``(first_frame, block)`` from an acquisition as it is written.

    Blocks are ``chunk`` frames, except that whatever has accumulated is handed
    over as soon as no new frames arrive: waiting for a full block would hold
    the last seconds of a sparse acquisition back indefinitely.

    ``stop_event`` is a :class:`threading.Event` that ends the stream early --
    what a viewer sets when its window is closed.  ``on_wait(seconds)`` is
    called on each idle poll with the time since the last new frame, for a
    progress line.
    """
    settings = settings or WatchSettings()
    first = _wait_for_file(path, settings, stop_event)
    if first is None:
        return

    index = 0                  # frames seen, across the whole series
    buffer: List[np.ndarray] = []
    buffer_start = start
    file_i, page_i = 0, 0
    last_new = time.monotonic()
    done = False
    # the very last page is held back while frames are still coming, but once
    # the acquisition has clearly ended nothing is going to follow it, so the
    # timeout triggers one more pass that does read it
    final = False

    while not _stopped(stop_event):
        files = _series_files(first)
        arrived = 0

        while file_i < len(files) and not done:
            is_last = file_i == len(files) - 1
            hold = settings.hold_last and is_last and not final
            images, page_i, exhausted = _read_new_pages(files[file_i], page_i,
                                                        hold_last=hold)
            for image in images:
                if (stop is not None and index >= stop) or _stopped(stop_event):
                    done = True
                    break
                if index >= start:
                    buffer.append(image)
                    if len(buffer) == chunk:
                        yield buffer_start, np.stack(buffer)
                        buffer_start += len(buffer)
                        buffer = []
                index += 1
                arrived += 1
            if done or not exhausted or is_last:
                break
            file_i += 1        # this file is complete and read; the next one
            page_i = 0

        if arrived and not (final or done):
            last_new = time.monotonic()
            continue           # there may be more waiting; look again at once

        if buffer:             # hand over a partial block rather than sit on it
            yield buffer_start, np.stack(buffer)
            buffer_start += len(buffer)
            buffer = []
        if final or done:
            break

        waited = time.monotonic() - last_new
        if waited >= settings.timeout:
            final = True
            continue
        if on_wait is not None:
            on_wait(waited)
        _sleep(settings.poll, stop_event)


def open_growing_stack(path, settings: Optional[WatchSettings] = None,
                       stop_event=None) -> ImageSource:
    """Open an acquisition as soon as its first frame is readable.

    The returned source describes the camera and the frame shape, which is what
    the pipeline needs to be set up; its ``n_frames`` is only what had been
    written by the time it was opened.  Use `watch_stack` for the frames.
    """
    settings = settings or WatchSettings()
    first = _wait_for_file(path, settings, stop_event)
    if first is None:
        raise TimeoutError(f"no readable TIFF appeared at {path} within "
                           f"{settings.appear_timeout} s")
    deadline = time.monotonic() + settings.appear_timeout
    while True:
        try:
            source = open_stack(first)
            if source.n_frames:
                return source
        except Exception:
            pass               # still being created; see the module docstring
        if time.monotonic() > deadline or _stopped(stop_event):
            raise TimeoutError(f"{first} has no readable frame yet")
        _sleep(settings.poll, stop_event)


# ------------------------------------------------------------------- internals
def _read_new_pages(path: Path, page_i: int, hold_last: bool
                    ) -> Tuple[List[np.ndarray], int, bool]:
    """Pages of ``path`` from ``page_i`` on, and whether the file is exhausted.

    The file is opened afresh every time: `tifffile` reads the page list once,
    so a handle kept open would never see the pages written after it.
    """
    try:
        with _quiet_tifffile(), tifffile.TiffFile(path) as tf:
            available = len(tf.pages)
            if hold_last:
                available -= 1
            images = []
            for i in range(page_i, max(available, page_i)):
                images.append(tf.pages[i].asarray())
            return images, page_i + len(images), page_i + len(images) >= available
    except Exception:
        # a file caught mid-write: nothing new this time round
        return [], page_i, False


def _wait_for_file(path, settings: WatchSettings, stop_event) -> Optional[Path]:
    """The base file of the series, once one exists."""
    path = Path(path)
    deadline = time.monotonic() + settings.appear_timeout
    while not _stopped(stop_event):
        if path.is_dir():
            found = sorted(path.glob("*.tif")) + sorted(path.glob("*.tiff"))
            if found:
                return _first_of_series(found)
        elif path.exists():
            return path
        if time.monotonic() > deadline:
            return None
        _sleep(settings.poll, stop_event)
    return None


@contextlib.contextmanager
def _quiet_tifffile():
    """Silence tifffile while reading a file that is being written.

    A half-written page chain makes it log an error per read; here that is the
    normal case and the answer is to look again, so the log would say nothing
    but "still writing", thousands of times over an acquisition.
    """
    logger = logging.getLogger("tifffile")
    level = logger.level
    logger.setLevel(logging.CRITICAL)
    try:
        yield
    finally:
        logger.setLevel(level)


def _stopped(stop_event) -> bool:
    return stop_event is not None and stop_event.is_set()


def _sleep(seconds: float, stop_event) -> None:
    """Sleep, but wake at once if the stream has been asked to stop."""
    if stop_event is None:
        time.sleep(seconds)
    else:
        stop_event.wait(seconds)
