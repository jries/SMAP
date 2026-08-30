"""Putting the stages together.

:class:`LocalizationEngine` is the actual driver: images are pushed in, ROIs
accumulate until there are enough to fit efficiently, and localizations come
out.  It never asks how many frames there will be, so the same object serves an
offline file and a live acquisition; :func:`fit_stack` is just a loop over a
file wrapped around it.

All state that has to survive between blocks -- the ROI buffer -- lives here,
so the stages themselves stay stateless functions.
"""

from __future__ import annotations

import queue
import threading
import time
from dataclasses import dataclass, field
from typing import Callable, Iterable, Iterator, List, Optional, Tuple

import numpy as np

from .camera import to_photons
from .detect import Candidates, PeakFinder
from .locs import Localizations, fit_to_localizations, to_nm, valid
from .metadata import CameraMetadata
from .psf import PSFModel
from .roi import ROIStack, cut_rois


@dataclass
class FitSettings:
    """Everything about the fit that is not the camera or the PSF model."""

    roisize: int = 13
    iterations: int = 50
    max_block_rois: int = 15000  # ROIs accumulated before a fit is run
    n_threads: int = 0  # 0 = one per core, for detection and fitting alike
    max_fit_distance: Optional[float] = None  # reject runaway fits, in pixels
    output_unit: str = "pixel"  # "pixel", "nm", or "pixel+nm"

    def block_rois(self) -> int:
        """SMAP scales the block with the ROI area; keep that behaviour."""
        return max(1, int(self.max_block_rois * 100 / self.roisize ** 2))


@dataclass
class LocalizationEngine:
    """Turns blocks of raw frames into localizations.

    ``push`` returns localizations whenever enough ROIs have accumulated,
    otherwise ``None``.  ``flush`` fits whatever is left.  For online use, call
    ``flush`` on a timer as well, so sparse data still yields results promptly.
    """

    camera: CameraMetadata
    finder: PeakFinder
    model: PSFModel
    settings: FitSettings = field(default_factory=FitSettings)
    mirror: Optional[bool] = None  # default: ask the PSF model

    _rois: List[np.ndarray] = field(default_factory=list, init=False)
    _candidates: List[Candidates] = field(default_factory=list, init=False)
    _buffered: int = field(default=0, init=False)
    stats: dict = field(default_factory=dict, init=False)

    def __post_init__(self):
        if self.mirror is None:
            self.mirror = bool(getattr(self.model, "mirror", False))
        self.stats = {"frames": 0, "candidates": 0, "rois": 0, "localizations": 0,
                      "dropped_at_border": 0, "rejected": 0,
                      "detect_seconds": 0.0, "fit_seconds": 0.0}

    # --------------------------------------------------------------- driving it
    def push(self, frames: np.ndarray, first_frame: int = 0
             ) -> Optional[Localizations]:
        """Process a block of raw frames ``(n, y, x)``."""
        t0 = time.perf_counter()
        photons = to_photons(frames, self.camera) / self.camera.excess_noise
        candidates, _ = self.finder(photons, first_frame=first_frame,
                                    n_threads=self.settings.n_threads)
        rois = cut_rois(photons, candidates, self.settings.roisize,
                        first_frame=first_frame, mirror=self.mirror)
        self.stats["detect_seconds"] += time.perf_counter() - t0
        self.stats["frames"] += len(photons)
        self.stats["candidates"] += len(candidates)
        self.stats["dropped_at_border"] += len(candidates) - len(rois)

        if len(rois):
            self._rois.append(rois.images)
            self._candidates.append(rois.candidates)
            self._buffered += len(rois)

        if self._buffered >= self.settings.block_rois():
            return self.flush()
        return None

    def flush(self) -> Optional[Localizations]:
        """Fit everything buffered so far."""
        if not self._rois:
            return None

        images = np.concatenate(self._rois) if len(self._rois) > 1 else self._rois[0]
        candidates = _concat_candidates(self._candidates)
        self._rois, self._candidates, self._buffered = [], [], 0

        stack = ROIStack(images, candidates, self.settings.roisize, self.mirror)
        t0 = time.perf_counter()
        result = self.model.fit(stack.images, iterations=self.settings.iterations,
                                n_threads=self.settings.n_threads)
        self.stats["fit_seconds"] += time.perf_counter() - t0
        self.stats["rois"] += len(stack)

        locs = fit_to_localizations(result, stack, self.model, self.camera)
        keep = valid(locs, self.settings.max_fit_distance)
        self.stats["rejected"] += int((~keep).sum())
        locs = locs[keep]
        self.stats["localizations"] += len(locs)

        if self.settings.output_unit != "pixel":
            self.camera.require("pixelsize_um")
            locs = to_nm(locs, float(self.camera.pixelsize_um) * 1000.0,
                         keep_pixels=self.settings.output_unit == "pixel+nm")
        return locs

    def __str__(self) -> str:
        s = self.stats
        return (f"{s['localizations']} localizations from {s['frames']} frames "
                f"({s['detect_seconds']:.1f} s detection, "
                f"{s['fit_seconds']:.1f} s fitting)")


def prefetch(source: Iterable, depth: int = 2) -> Iterator:
    """Yield from ``source`` while a background thread reads ahead.

    Reading frames from disk and processing them are independent, and both
    release the GIL for most of their work, so overlapping them costs one
    thread and hides the read time almost entirely.  ``depth`` bounds how many
    blocks may be held in memory at once, which keeps this usable for a long
    acquisition.
    """
    if depth <= 0:
        yield from source
        return

    items: "queue.Queue" = queue.Queue(maxsize=depth)
    sentinel = object()

    def reader():
        try:
            for item in source:
                items.put(item)
        except BaseException as error:  # re-raised in the consumer
            items.put(error)
        finally:
            items.put(sentinel)

    thread = threading.Thread(target=reader, daemon=True)
    thread.start()
    try:
        while True:
            item = items.get()
            if item is sentinel:
                break
            if isinstance(item, BaseException):
                raise item
            yield item
    finally:
        thread.join(timeout=1.0)


def _concat_candidates(parts: List[Candidates]) -> Candidates:
    if len(parts) == 1:
        return parts[0]
    return Candidates(
        frame=np.concatenate([p.frame for p in parts]),
        x=np.concatenate([p.x for p in parts]),
        y=np.concatenate([p.y for p in parts]),
        value=np.concatenate([p.value for p in parts]),
    )


def fit_stack(frames: Iterable[Tuple[int, np.ndarray]], camera: CameraMetadata,
              finder: PeakFinder, model: PSFModel,
              settings: Optional[FitSettings] = None,
              sink: Optional[Callable[[Localizations], None]] = None,
              progress: Optional[Callable[[LocalizationEngine], None]] = None,
              read_ahead: int = 2
              ) -> Tuple[Localizations, LocalizationEngine]:
    """Run the pipeline over an iterable of ``(first_frame, block)``.

    ``sink`` is called with each finished block of localizations -- pass the
    ``append`` method of a :class:`~smapfit.io.hdf5.LocalizationWriter` to
    stream to disk instead of keeping everything in memory.

    ``read_ahead`` blocks are fetched from ``frames`` by a background thread;
    set it to 0 to read strictly in step with processing.
    """
    engine = LocalizationEngine(camera, finder, model, settings or FitSettings())
    collected = Localizations()

    def emit(locs: Optional[Localizations]) -> None:
        if locs is None or len(locs) == 0:
            return
        if sink is not None:
            sink(locs)
        else:
            collected.extend(locs)

    for first_frame, block in prefetch(frames, read_ahead):
        emit(engine.push(block, first_frame))
        if progress is not None:
            progress(engine)
    emit(engine.flush())

    return collected, engine


def provenance(camera: CameraMetadata, finder: PeakFinder, model: PSFModel,
               settings: FitSettings, source=None) -> dict:
    """A record of how a result was produced, for the output file."""
    info = {
        "source": str(source) if source is not None else None,
        "camera": camera.to_dict(),
        "detection": {"filter": str(finder.filter), "cutoff": str(finder.cutoff)},
        "psf_model": str(model),
        "fit": {"roisize": settings.roisize, "iterations": settings.iterations,
                "rois_per_block": settings.block_rois(),
                "max_fit_distance": settings.max_fit_distance,
                "output_unit": settings.output_unit},
        "smapfit_version": _version(),
    }
    calibration = getattr(model, "calibration", None)
    if calibration is not None:
        info["calibration"] = {
            "file": str(calibration.source),
            "dz_nm": calibration.dz, "z0": calibration.z0,
            "grid": list(calibration.shape),
            "mirrored": bool(calibration.em_mirror),
            "em_on": calibration.em_on,
        }
    return info


def _version() -> str:
    try:
        from importlib.metadata import version
        return version("smapfit")
    except Exception:
        return "0.1.0dev"
