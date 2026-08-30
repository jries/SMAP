"""Camera information from SMAP's ``settings/*_cameras.mat``.

Micro-Manager does not record the e-/ADU conversion, which depends on the
readout mode.  SMAP keeps that in a settings file: each camera has a list of
*states* (readout port, readout rate, gain), and each state carries the
conversion for that mode.  The conversion is the only measured **value** this
module takes from the file.

The file also describes, per camera, *which metadata key* means what -- for an
Evolve the EM mode is ``Evolve512-Port != 'Normal'``, for an Andor it is
``Andor-Output_Amplifier != 'Conventional'``.  Those are interpretation rules,
not values, and we use them so that EM settings are read correctly from the
image metadata of any camera in the file.

This is a stopgap for compatibility with the existing lab settings file; the
intended long-term source is a plain config file.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import scipy.io

# metadata-derived fields we are willing to interpret via the preset rules
_RULE_FIELDS = {"EMon", "emgain", "exposure", "offset"}


@dataclass
class MetadataRule:
    """How to read one quantity out of the Micro-Manager metadata."""

    key: str  # the metadata key to read
    expression: str  # a small MATLAB expression applied to its value

    def apply(self, mm: Dict[str, object]):
        if self.key not in mm:
            return None
        return _eval_expr(self.expression, str(mm[self.key]))


def _eval_expr(expr: str, x: str):
    """Evaluate the handful of MATLAB expressions used in the settings file."""
    expr = (expr or "").strip()
    if not expr or expr == "[]":
        return x
    if expr in ("str2double(X)", "str2num(X)"):
        return _as_float(x)
    if expr == "str2double(X)>0":
        v = _as_float(x)
        return None if v is None else v > 0
    m = re.fullmatch(r"\(?(~?)strcmp\(X,'([^']*)'\)\)?", expr)
    if m:
        equal = x.strip() == m.group(2)
        return (not equal) if m.group(1) else equal
    return x


@dataclass
class CameraState:
    """One readout mode: the metadata that identifies it, and its conversion."""

    match: Dict[str, str]
    conversion: Optional[float]
    offset: Optional[float] = None

    def matches(self, mm: Dict[str, object]) -> bool:
        conditions = {k: v for k, v in self.match.items()
                      if k and k != "select" and v != ""}
        if not conditions:
            return False
        return all(str(mm.get(k, "\0")).strip() == v.strip()
                   for k, v in conditions.items())


@dataclass
class CameraEntry:
    name: str
    id_key: str
    id_value: str
    states: List[CameraState]
    rules: Dict[str, MetadataRule] = field(default_factory=dict)

    def identifies(self, mm: Dict[str, object]) -> bool:
        if not self.id_key:
            return False
        value = mm.get(self.id_key)
        return value is not None and str(value).strip() == self.id_value.strip()


def _s(v) -> str:
    """MATLAB cell entry -> plain string ('' for empty)."""
    if isinstance(v, np.ndarray):
        return "" if v.size == 0 else str(v.reshape(-1)[0])
    return "" if v is None else str(v)


def _as_float(v) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


class CameraPresets:
    """Cameras defined in a SMAP ``*_cameras.mat`` file."""

    def __init__(self, cameras: List[CameraEntry], source: Optional[Path] = None):
        self.cameras = cameras
        self.source = source

    @classmethod
    def load(cls, path) -> "CameraPresets":
        path = Path(path)
        mat = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
        camtab = np.atleast_2d(mat["camtab"])
        raw = np.atleast_1d(mat["cameras"])

        cameras = []
        for row, cam in zip(camtab, raw):
            states = []
            for st in np.atleast_1d(getattr(cam, "state", [])):
                if not hasattr(st, "_fieldnames"):
                    continue
                values = {_s(r[0]): _s(r[1])
                          for r in np.atleast_2d(getattr(st, "par", []))}
                states.append(CameraState(
                    match={_s(r[0]): _s(r[1])
                           for r in np.atleast_2d(getattr(st, "defpar", []))},
                    conversion=_as_float(values.get("conversion")),
                    offset=_as_float(values.get("offset")),
                ))
            rules = {}
            for r in np.atleast_2d(cam.par):
                name, mode = _s(r[0]), _s(r[1])
                if mode == "metadata" and name in _RULE_FIELDS:
                    rules[name] = MetadataRule(key=_s(r[3]), expression=_s(r[5]))
            cameras.append(CameraEntry(name=_s(row[0]), id_key=_s(row[1]),
                                       id_value=_s(row[2]), states=states,
                                       rules=rules))
        return cls(cameras, path)

    def identify(self, mm: Dict[str, object]) -> Optional[CameraEntry]:
        for cam in self.cameras:
            if cam.identifies(mm):
                return cam
        return None

    def state_for(self, mm: Dict[str, object]) -> Optional[CameraState]:
        cam = self.identify(mm)
        if cam is None:
            return None
        for state in cam.states:
            if state.matches(mm):
                return state
        return None

    def conversion_for(self, mm: Dict[str, object]) -> Optional[float]:
        """e-/ADU for the readout mode described by ``mm``, if it can be found."""
        state = self.state_for(mm)
        return None if state is None else state.conversion

    def interpret(self, mm: Dict[str, object]) -> dict:
        """Read camera settings out of image metadata using this camera's rules.

        Returns a dict that may contain ``em_on``, ``emgain``, ``exposure_ms``
        and ``offset`` -- all read from the *image metadata*, using the keys and
        expressions the settings file defines for this camera -- plus
        ``conversion`` and ``state_offset``, which come from the settings file
        itself.  Absent entries mean the metadata did not provide them.
        """
        cam = self.identify(mm)
        if cam is None:
            return {}

        out: dict = {"camera_preset": cam.name}
        mapping = {"EMon": "em_on", "emgain": "emgain",
                   "exposure": "exposure_ms", "offset": "offset"}
        for name, rule in cam.rules.items():
            value = rule.apply(mm)
            if value is not None:
                out[mapping[name]] = value

        if isinstance(out.get("em_on"), str):
            out["em_on"] = out["em_on"].strip().lower() == "true"

        state = self.state_for(mm)
        if state is not None:
            out["conversion"] = state.conversion
            out["state_offset"] = state.offset
        return out

    def describe(self, mm: Dict[str, object]) -> str:
        cam = self.identify(mm)
        if cam is None:
            return "camera not found in presets"
        conv = self.conversion_for(mm)
        if conv is None:
            return f"{cam.name}: no matching readout state"
        return f"{cam.name}: conversion {conv} e-/ADU"
