"""Writing and reading localization tables as HDF5.

One file holds everything: the localization columns under ``/locs`` and the
full provenance -- camera metadata, every stage's settings, the calibration
used -- as JSON in the file attributes.

The writer appends block by block and flushes as it goes, so results are on
disk and readable while an acquisition is still running, and an interrupted run
keeps everything up to the last block.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional

import h5py
import numpy as np

from ..locs import Localizations

FORMAT_VERSION = 1


class LocalizationWriter:
    """Append-only HDF5 writer for localizations.

    Usage::

        with LocalizationWriter("out.h5", metadata=...) as writer:
            writer.append(locs)
    """

    def __init__(self, path, metadata: Optional[Dict[str, object]] = None,
                 chunk: int = 8192, overwrite: bool = True,
                 compression: Optional[str] = "lzf"):
        self.path = Path(path)
        self.chunk = int(chunk)
        # lzf compresses about five times faster than gzip for ~15% more space,
        # which matters because this runs while frames are still being fitted
        self.compression = compression
        mode = "w" if overwrite else "w-"
        self._file = h5py.File(self.path, mode)
        self._group = self._file.create_group("locs")
        self._datasets: Dict[str, h5py.Dataset] = {}
        self._metadata: Dict[str, object] = {}
        self._n = 0
        self._file.attrs["format"] = "smappy-localizations"
        self._file.attrs["format_version"] = FORMAT_VERSION
        if metadata:
            self.set_metadata(metadata)

    # ------------------------------------------------------------------ writing
    def set_metadata(self, metadata: Dict[str, object]) -> None:
        """Store provenance as JSON; may be called again to add more."""
        self._metadata.update(metadata)
        self._write_metadata()

    def _write_metadata(self) -> None:
        self._file.attrs["metadata"] = json.dumps(self._metadata,
                                                  default=_json_default, indent=1)

    def append(self, locs: Localizations) -> None:
        if len(locs) == 0:
            return
        if not self._datasets:
            self._create(locs)
            # the table itself knows its units; keep that with the provenance
            if locs.metadata:
                self.set_metadata(locs.metadata)
        elif set(self._datasets) != set(locs.keys()):
            raise ValueError(
                "columns changed between blocks: "
                f"{sorted(set(self._datasets) ^ set(locs.keys()))}")

        n = len(locs)
        for name, dataset in self._datasets.items():
            dataset.resize(self._n + n, axis=0)
            dataset[self._n:self._n + n] = locs[name]
        self._n += n
        self._file.attrs["n_localizations"] = self._n
        self._file.flush()

    def _create(self, locs: Localizations) -> None:
        for name in sorted(locs.keys()):
            column = np.asarray(locs[name])
            self._datasets[name] = self._group.create_dataset(
                name, shape=(0,), maxshape=(None,), dtype=column.dtype,
                chunks=(self.chunk,), compression=self.compression)

    # ------------------------------------------------------------------ closing
    def close(self) -> None:
        if self._file:
            self._file.attrs["n_localizations"] = self._n
            self._file.close()
            self._file = None

    def __enter__(self) -> "LocalizationWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __len__(self) -> int:
        return self._n


def save_localizations(path, locs: Localizations,
                       metadata: Optional[Dict[str, object]] = None) -> Path:
    """Write a complete table in one go."""
    with LocalizationWriter(path, metadata or locs.metadata) as writer:
        writer.append(locs)
    return Path(path)


def load_localizations(path) -> Localizations:
    """Read a table written by :class:`LocalizationWriter`."""
    with h5py.File(path, "r") as f:
        columns = {name: f["locs"][name][()] for name in f["locs"]}
        metadata = {}
        if "metadata" in f.attrs:
            metadata = json.loads(f.attrs["metadata"])
    return Localizations(columns, metadata)


def _json_default(obj):
    """Make numpy values and paths JSON-serializable."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)
