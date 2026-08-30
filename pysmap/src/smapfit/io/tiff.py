"""Reading Micro-Manager / OME TIFF stacks and their metadata.

Micro-Manager writes one acquisition as a series of files
(``..._MMStack_Default.ome.tif``, ``..._1.ome.tif``, ...) and records in the
summary metadata how many frames were *planned*.  Acquisitions are usually
stopped early, so the declared frame count is larger than what was written.
We therefore never trust the declared length: frames are taken from the pages
actually present in each file, in acquisition order.

The loader yields chunks of frames so that downstream stages can work on
batches, and never requires the whole stack to be in memory or the source to
be finite -- the same interface works for online analysis later.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import tifffile

from ..metadata import CameraMetadata


@dataclass
class ImageSource:
    """A Micro-Manager acquisition, possibly split over several files."""

    files: List[Path]
    shape: Tuple[int, int]  # (height, width) of one frame
    dtype: np.dtype
    n_frames: int  # frames actually written
    n_frames_declared: Optional[int]  # what MM planned, if it said
    mm_metadata: Dict[str, object]  # per-plane metadata of the first frame
    summary: Dict[str, object]

    def __len__(self) -> int:
        return self.n_frames

    # ------------------------------------------------------------------ frames
    def frames(self, chunk: int = 100, start: int = 0,
               stop: Optional[int] = None) -> Iterator[Tuple[int, np.ndarray]]:
        """Yield ``(first_frame_index, block)`` with block shape (n, y, x).

        Frame indices are zero-based and continuous across the files of the
        acquisition.
        """
        index = 0
        buffer: List[np.ndarray] = []
        buffer_start = start

        for path in self.files:
            with tifffile.TiffFile(path) as tf:
                for page in tf.pages:
                    if stop is not None and index >= stop:
                        break
                    if index < start:
                        index += 1
                        continue
                    buffer.append(page.asarray())
                    if len(buffer) == chunk:
                        yield buffer_start, np.stack(buffer)
                        buffer_start += len(buffer)
                        buffer = []
                    index += 1
            if stop is not None and index >= stop:
                break

        if buffer:
            yield buffer_start, np.stack(buffer)

    def frame(self, index: int) -> np.ndarray:
        """Read a single frame (convenience for previews and tests)."""
        for start, block in self.frames(chunk=1, start=index, stop=index + 1):
            return block[0]
        raise IndexError(f"frame {index} beyond end of stack ({self.n_frames})")


def open_stack(path) -> ImageSource:
    """Open a Micro-Manager TIFF acquisition, including its continuation files."""
    path = Path(path)
    if path.is_dir():
        candidates = sorted(path.glob("*.tif")) + sorted(path.glob("*.tiff"))
        if not candidates:
            raise FileNotFoundError(f"no TIFF files in {path}")
        path = _first_of_series(candidates)

    files = _series_files(path)

    n_frames = 0
    for f in files:
        with tifffile.TiffFile(f) as tf:
            n_frames += len(tf.pages)

    with tifffile.TiffFile(files[0]) as tf:
        page = tf.pages[0]
        shape, dtype = tuple(page.shape), page.dtype
        summary = {}
        mm = tf.micromanager_metadata or {}
        if isinstance(mm.get("Summary"), dict):
            summary = mm["Summary"]
        plane = {}
        tag = page.tags.get("MicroManagerMetadata")
        if tag is not None and isinstance(tag.value, dict):
            plane = tag.value

    declared = summary.get("Frames")
    return ImageSource(
        files=files, shape=shape, dtype=dtype, n_frames=n_frames,
        n_frames_declared=int(declared) if declared else None,
        mm_metadata=plane, summary=summary,
    )


def _first_of_series(paths: List[Path]) -> Path:
    """Pick the base file of an MM series (the one without a ``_<n>`` suffix)."""
    base = [p for p in paths if not re.search(r"_\d+\.ome\.tiff?$", p.name)]
    return (base or paths)[0]


def _series_files(first: Path) -> List[Path]:
    """All files of a Micro-Manager series, in acquisition order."""
    name = first.name
    m = re.match(r"(?P<stem>.+?)(?:_(?P<n>\d+))?\.ome\.tiff?$", name)
    if not m:
        return [first]

    stem = m.group("stem")
    found = {}
    for candidate in first.parent.glob(f"{stem}*.ome.tif*"):
        mm = re.match(rf"{re.escape(stem)}(?:_(\d+))?\.ome\.tiff?$", candidate.name)
        if mm:
            found[int(mm.group(1) or 0)] = candidate
    return [found[k] for k in sorted(found)] or [first]


# --------------------------------------------------------------------- metadata
def metadata_from_stack(source: ImageSource,
                        presets=None) -> CameraMetadata:
    """Extract what the image file knows about the camera.

    ``presets`` is an optional :class:`~smapfit.io.cameras_mat.CameraPresets`.
    It supplies the e-/ADU conversion, which Micro-Manager does not record, and
    the per-camera rules that say which metadata key carries the EM mode, the
    EM gain and the offset -- those differ between an Evolve and an iXon.

    The offset is taken from the image metadata.  Some cameras (e.g. the iXon)
    do not report it at all; in that case the value stored in the settings file
    for the matching readout mode is used and a warning is issued.  Values that
    remain unknown stay ``None`` and must be supplied by the user.
    """
    mm = source.mm_metadata
    device = str(mm.get("Core-Camera") or mm.get("Camera") or "").strip()

    def dev(key, default=None):
        return mm.get(f"{device}-{key}", default) if device else default

    info = presets.interpret(mm) if presets is not None else {}

    # EM: use the camera's own rule when we have one, else a generic guess
    if "em_on" in info:
        em_on = bool(info["em_on"])
    else:
        port = dev("Port") or dev("Output_Amplifier")
        em_on = bool(port is not None
                     and str(port).strip() not in ("Normal", "Conventional"))
    emgain = _to_float(info.get("emgain")) or _to_float(dev("MultiplierGain")) \
        or _to_float(dev("Gain")) or 1.0

    offset = _to_float(info.get("offset"))
    if offset is None:
        offset = _to_float(dev("Offset"))
    if offset is None and info.get("state_offset") is not None:
        offset = _to_float(info["state_offset"])
        warnings.warn(
            f"camera '{device or '?'}' does not report an offset in the image "
            f"metadata; using {offset} ADU from the settings file "
            f"({info.get('camera_preset')}). Set it explicitly to be sure.",
            stacklevel=2)

    pixelsize = _to_float(mm.get("PixelSizeUm"))
    if not pixelsize:  # MM reports 0.0 when the pixel size is not calibrated
        pixelsize = None

    return CameraMetadata(
        conversion=_to_float(info.get("conversion")),
        offset=offset,
        pixelsize_um=pixelsize,
        em_on=em_on,
        emgain=emgain if em_on else 1.0,
        roi=_parse_roi(mm.get("ROI")),
        exposure_ms=(_to_float(info.get("exposure_ms"))
                     or _to_float(mm.get("Exposure-ms"))
                     or _to_float(dev("Exposure"))),
        camera_name=device or None,
    )


def camera_metadata(source: ImageSource, presets=None, overrides=None,
                    require: bool = True) -> CameraMetadata:
    """Camera metadata for a stack: file metadata, then user overrides.

    ``overrides`` may be a :class:`~smapfit.metadata.CameraMetadata`, a dict, or
    a path to a YAML file.  Anything set there wins over the image metadata.
    With ``require=True`` the result is checked for completeness, so a missing
    pixel size fails here rather than silently producing wrong nm coordinates.
    """
    meta = metadata_from_stack(source, presets)

    if overrides is not None:
        if isinstance(overrides, CameraMetadata):
            user = overrides
        elif isinstance(overrides, dict):
            user = CameraMetadata.from_dict(overrides)
        else:
            user = CameraMetadata.from_yaml(overrides)
        meta = meta.merged_with(user)

    if require:
        meta.require()
    return meta


def _to_float(v) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _parse_roi(value) -> Optional[Tuple[int, int, int, int]]:
    """Micro-Manager writes the ROI as ``"x-y-width-height"``.

    The separator is a hyphen, so the numbers must not be read as signed.
    """
    if value is None:
        return None
    parts = re.findall(r"\d+", str(value))
    if len(parts) != 4:
        return None
    return tuple(int(p) for p in parts)
