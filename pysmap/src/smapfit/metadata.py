"""Camera metadata needed to turn raw camera counts into photons.

Only what the fitting pipeline actually needs is kept.  Values can come from
three places, in increasing priority:

1. the image file itself (Micro-Manager / OME tags)
2. a camera preset (:mod:`smapfit.io.cameras_mat`, SMAP's ``*_cameras.mat``)
3. explicit user input -- a dict or a YAML file

Missing required values raise a clear error rather than silently defaulting.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
from typing import Optional, Sequence

REQUIRED = ("conversion", "offset", "pixelsize_um")


@dataclass
class CameraMetadata:
    """Camera parameters for ADU -> photon conversion and nm coordinates."""

    conversion: Optional[float] = None  # e- per ADU
    offset: Optional[float] = None  # camera baseline, ADU
    pixelsize_um: Optional[float] = None  # effective pixel size in the sample
    em_on: Optional[bool] = None  # EM gain used (EMCCD); None = unspecified
    emgain: Optional[float] = None  # EM multiplication gain
    roi: Optional[Sequence[int]] = None  # (x, y, width, height) on the chip
    exposure_ms: Optional[float] = None
    camera_name: Optional[str] = None
    comment: str = ""

    # ---------------------------------------------------------------- factors
    @property
    def adu_to_photons(self) -> float:
        """Multiplicative factor applied to (counts - offset)."""
        self.require("conversion")
        if self.is_em:
            if not self.emgain:
                raise ValueError("em_on is set but emgain is zero or unset")
            return float(self.conversion) / float(self.emgain)
        return float(self.conversion)

    @property
    def is_em(self) -> bool:
        """Whether EM amplification was used (unspecified counts as off)."""
        return bool(self.em_on)

    @property
    def excess_noise(self) -> float:
        """EMCCD excess-noise factor.

        The fitter assumes pure Poisson noise.  EM amplification roughly
        doubles the variance, which is handled outside the fitter by dividing
        the photon counts by this factor before fitting and multiplying
        photons and background by it afterwards.
        """
        return 2.0 if self.is_em else 1.0

    @property
    def roi_offset(self) -> tuple:
        """(x, y) chip offset of the image, used for absolute nm coordinates."""
        if self.roi is None:
            return (0, 0)
        return (int(self.roi[0]), int(self.roi[1]))

    # ---------------------------------------------------------------- merging
    def merged_with(self, other: "CameraMetadata") -> "CameraMetadata":
        """Return a copy where set values of ``other`` override ours."""
        upd = {f.name: getattr(other, f.name) for f in fields(other)
               if not _is_unset(getattr(other, f.name), f.name)}
        return replace(self, **upd)

    def require(self, *names: str) -> None:
        missing = [n for n in (names or REQUIRED) if getattr(self, n) is None]
        if missing:
            raise ValueError(
                "missing camera metadata: " + ", ".join(missing) +
                ". Supply it via a camera preset or explicitly, e.g. "
                "CameraMetadata(conversion=6.7, offset=398.6, pixelsize_um=0.127)"
            )

    def to_dict(self) -> dict:
        d = asdict(self)
        d["roi"] = None if self.roi is None else [int(v) for v in self.roi]
        return d

    # ------------------------------------------------------------ contructors
    @classmethod
    def from_dict(cls, d: dict) -> "CameraMetadata":
        known = {f.name for f in fields(cls)}
        unknown = set(d) - known
        if unknown:
            raise ValueError(f"unknown camera metadata fields: {sorted(unknown)}")
        return cls(**d)

    @classmethod
    def from_yaml(cls, path) -> "CameraMetadata":
        import yaml
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        if "camera" in data:  # allow a nested section in a larger config
            data = data["camera"]
        return cls.from_dict(data)

    def to_yaml(self, path) -> None:
        import yaml
        with open(path, "w") as fh:
            yaml.safe_dump(self.to_dict(), fh, sort_keys=False)

    def __str__(self) -> str:
        px = "?" if self.pixelsize_um is None else f"{self.pixelsize_um:g} um"
        em = (f"EM on, gain {self.emgain:g}" if self.is_em
              else ("EM off" if self.em_on is not None else "EM unspecified"))
        return (f"{self.camera_name or 'camera'}: conversion="
                f"{self.conversion} e-/ADU, offset={self.offset} ADU, "
                f"pixel {px}, {em}, roi={list(self.roi) if self.roi is not None else None}")


def _is_unset(value, name: str) -> bool:
    """Values that mean 'not specified' when merging.

    Every optional field defaults to ``None``, so a user config that omits
    ``em_on`` cannot accidentally switch EM off.
    """
    if value is None:
        return True
    if name == "comment" and value == "":
        return True
    return False


def load_camera_metadata(path) -> CameraMetadata:
    """Load user-supplied metadata from a YAML file."""
    return CameraMetadata.from_yaml(Path(path))
