"""Drift correction: estimate with COMET, apply to every localization.

Sample drift is a property of the *acquisition*, not of a subset of it, so the
two halves of drift correction want different data:

* **Estimating** it works best on a clean subset -- bright, well-fitted
  localizations -- because a badly localized point contributes noise to the
  overlap cost and nothing else.  That is exactly what the display filter
  already selects, so `estimate_drift` takes the same `select` a render does.
* **Applying** it is a coordinate correction and must reach *all*
  localizations, including the ones the filter hides.  A filter is a view;
  changing it later must not require re-running the correction.

The estimator is COMET (Cost-function Optimized Maximal Overlap drift
EsTimation, https://github.com/gpufit/Comet), used as a library -- it maximises
the spatiotemporal overlap of localizations across time windows, which needs no
fiducials and no reference structure.  Its result is a per-frame drift table
whose row index *is* the frame number, so applying it is one fancy-index per
coordinate.

Not wrapped: COMET's own file I/O, plotting and molecule-set output.  We hand it
an array and get an array back; reading and writing stays with `io.hdf5`, so the
corrected file is an ordinary smapfit file that the viewer opens unchanged.

Units: COMET works in nm.  A table in pixels is converted for the estimate and
the drift is divided back by the pixel size when applied, so the corrected table
keeps whatever units it had.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Optional

import numpy as np

from .locs import Localizations

# the position columns a drift correction touches, and the drift axis each
# follows.  Deliberately not `peak_*`: those are detection positions, kept as
# the record of what was found where, not measurements of the structure.
_POSITION_COLUMNS = {"x_nm": 0, "x_pix": 0,
                     "y_nm": 1, "y_pix": 1,
                     "z_nm": 2}


@dataclass
class DriftSettings:
    """COMET parameters.  The defaults suit a typical SMLM movie.

    ``segmentation_mode`` sets what ``segmentation_var`` counts:
    0 = number of time windows, 1 = localizations per window,
    2 = frames per window (the default).

    A time window must hold enough localizations to define the structure --
    a few thousand is a good target -- and the movie must give at least a few
    tens of windows for the drift curve to have any time resolution.
    """

    segmentation_mode: int = 2
    segmentation_var: int = 500
    max_drift_nm: float = 300.0      # also the neighbour-search radius
    target_sigma_nm: float = 30.0    # refinement stops at this length scale
    initial_sigma_nm: Optional[float] = None   # default: max_drift_nm / 3
    boxcar_width: int = 1            # temporal smoothing between iterations
    interpolation: str = "cubic"     # "cubic" or "catmull-rom"
    max_locs_per_segment: Optional[int] = None   # cap, to bound memory
    # L-BFGS-B relative tolerance.  COMET's own is 1e3*eps (~2e-13), which buys
    # nothing here: loosening it to 1e-7 costs 0.14 nm rms (0.9 nm worst case)
    # on the drift curve -- far below any localization precision -- and cuts the
    # cost evaluations from 423 to 154.  Raise it towards 1e-13 to reproduce
    # upstream COMET exactly.
    optimizer_ftol: float = 1e-7
    # Compute the Gaussian from a table plus a polynomial and skip pairs beyond
    # 6 sigma: 1.3x faster, gradient accurate to 4e-4 relative.  See NOTES.
    approximate_kernel: bool = True
    # A looser tolerance for the coarse sigma steps.  None (the default) means
    # the same as `optimizer_ftol`: measured, loosening the coarse steps is a
    # bad trade -- 2.3x faster for tens of nm on individual windows.  See NOTES.
    optimizer_ftol_coarse: Optional[float] = None
    backend: Optional[str] = None    # "cuda", "torch", "cpu"; None = fastest
    # Estimate from *grouped* localizations: one entry per emitter blink rather
    # than one per frame.  20-35x faster (it removes exactly the dense, very
    # short-range pairs, which is where the pair count grows fastest) and, on
    # the clathrin dataset, recovers 96% of the improvement the full estimate
    # makes to the image.  Off by default: the full estimate is the accurate
    # one, and is now fast enough to be the thing you run.
    # Check each time window's estimate against the overlap it would have with
    # no drift correction at all, and discard the ones that do not beat it --
    # the spline then bridges them.  Costs one extra pass over the pairs.
    quality_control: bool = False   # time-window fit only, not the spline
    # A window is discarded when its fitted drift improves its own overlap by
    # less than this, relative to no correction at all.  0.0 keeps every window
    # whose estimate is better than nothing, which is the honest bar.
    min_lift: float = 0.0
    group: bool = True
    group_dx_nm: float = 50.0        # linking radius
    group_dt: int = 1                # frames a blink may be missing
    # Estimate in two passes: a grouped one to get within a few nm, then an
    # ungrouped one over a small radius to refine.  The pair count is what costs
    # time, and once the coarse pass has taken the drift out, a 30 nm radius
    # keeps a seventh of the pairs -- 18 s instead of 153 s for 99% of the
    # improvement.  The small radius also bounds the fine pass, which is what
    # keeps it from producing the runaway window the single pass has here.
    two_stage: bool = False
    two_stage_radius_nm: float = 30.0
    # Fit the drift as a cubic B-spline in time instead of one free vector per
    # time window.  Drift is smooth -- slow creep, occasionally a small jump --
    # so a curve with a coefficient every `spline_knot_frames` frames has far
    # fewer parameters than a window per 500 frames, and each coefficient is
    # constrained by every localization near it in time rather than by one
    # window's worth.  There is no time window and no interpolation afterwards:
    # every localization is placed at its own frame, and the fitted curve *is*
    # the per-frame drift.  B-splines are variation-diminishing, so the curve
    # cannot overshoot the way an interpolating cubic spline through noisy
    # window estimates does.
    spline: bool = True
    # 2000 frames was the knee on the clathrin dataset (26 coefficients over
    # 46 005 frames); finer knots are noisier without being sharper, and there
    # was no sign of over-smoothing even at 4000.  A dataset with faster drift
    # wants a smaller spacing, which shows up as out-of-sample sharpness getting
    # worse.  Quality control does not apply here -- there are no time windows.
    spline_knot_frames: int = 2000   # spacing of the coefficients, in frames
    spline_penalty: float = 0.0      # extra roughness penalty on 2nd differences
    use_z: Optional[bool] = None     # None: 3D if the table has z_nm


class Drift:
    """Per-frame drift in nm, indexed by frame number.

    ``drift[f]`` is ``(dx, dy, dz)``: how far the sample had moved in frame
    ``f``, so correcting a localization *subtracts* it.
    """

    def __init__(self, drift_nm: np.ndarray, settings: Optional[DriftSettings] = None,
                 n_used: Optional[int] = None,
                 flagged_windows: Optional[np.ndarray] = None):
        drift = np.asarray(drift_nm, dtype=np.float64)
        if drift.ndim != 2 or drift.shape[1] != 3:
            raise ValueError(f"drift must be (n_frames, 3), got {drift.shape}")
        self.drift = drift
        self.settings = settings
        self.n_used = n_used  # localizations the estimate was made from
        # time windows quality control discarded; the curve is interpolated
        # across them, so they are a caveat on the result, not a hole in it
        self.flagged_windows = flagged_windows

    def __len__(self) -> int:
        return self.drift.shape[0]

    @property
    def frames(self) -> np.ndarray:
        return np.arange(len(self))

    def __str__(self) -> str:
        span = self.drift.max(axis=0) - self.drift.min(axis=0)
        text = (f"drift over {len(self)} frames, range "
                f"x {span[0]:.0f} nm, y {span[1]:.0f} nm, z {span[2]:.0f} nm")
        if self.flagged_windows is not None and len(self.flagged_windows):
            text += (f"; {len(self.flagged_windows)} time window(s) flagged and "
                     f"interpolated over: {list(self.flagged_windows)}")
        return text

    def apply(self, locs: Localizations,
              pixelsize_nm: Optional[float] = None) -> Localizations:
        """Return a copy of the whole table with the drift subtracted."""
        frames = np.asarray(locs["frame"], dtype=np.int64)
        if len(frames) and (frames.min() < 0 or frames.max() >= len(self)):
            raise ValueError(
                f"the table spans frames {frames.min()}..{frames.max()} but the "
                f"drift covers 0..{len(self) - 1}; it was estimated from a "
                f"different acquisition")

        columns = dict(locs.columns)
        corrected = []
        for name, axis in _POSITION_COLUMNS.items():
            if name not in columns:
                continue
            shift = self.drift[frames, axis]
            if name.endswith("_pix"):
                shift = shift / _pixelsize_nm(locs, pixelsize_nm)
            values = np.asarray(columns[name], dtype=np.float64) - shift
            columns[name] = values.astype(np.asarray(columns[name]).dtype)
            corrected.append(name)

        metadata = dict(locs.metadata)
        metadata["drift_correction"] = {
            "method": "comet",
            "columns": corrected,
            "n_frames": len(self),
            "n_localizations_used": self.n_used,
            "flagged_windows": (None if self.flagged_windows is None
                                else [int(w) for w in self.flagged_windows]),
            "settings": asdict(self.settings) if self.settings else None,
        }
        return Localizations(columns, metadata)

    def plot(self, ax=None):
        """Drift vs frame, the standard sanity check."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        for axis, label in enumerate("xyz"):
            if np.any(self.drift[:, axis]):
                ax.plot(self.frames, self.drift[:, axis], label=label)
        ax.set_xlabel("frame")
        ax.set_ylabel("drift (nm)")
        ax.legend()
        return ax


# --------------------------------------------------------------------- running
def estimate_drift(locs: Localizations, settings: Optional[DriftSettings] = None,
                   select=None, pixelsize_nm: Optional[float] = None,
                   display: bool = False) -> Drift:
    """Estimate drift from the selected localizations.

    ``select`` is what to estimate *from*: a `LocFilter`, a boolean mask, an
    index array, or None for everything.  The result covers every frame of the
    input table, not only the selected ones.
    """
    from comet import best_backend, comet_run_kd  # optional; see NOTES.md
    from comet.core import drift_optimizer
    from comet.core.qc_utils import flag_flawed_segments_by_lift
    from comet.core.cpu_wrapper import (cpu_wrapper_chunked_approx,
                                        cpu_wrapper_chunked_fast)

    settings = settings or DriftSettings()
    if settings.quality_control and settings.spline:
        raise ValueError("quality control flags time windows, and the spline "
                         "fit has none; pass spline=False to use it")
    if settings.spline:
        return _estimate_spline(locs, settings, select, pixelsize_nm, display)
    if settings.two_stage:
        return _two_stage(locs, settings, select, pixelsize_nm, display)

    # Opt in to the two changes made to the vendored COMET (see NOTES.md); its
    # own defaults are untouched for anyone importing it directly.
    drift_optimizer.CPU_WRAPPER = (cpu_wrapper_chunked_approx
                                   if settings.approximate_kernel
                                   else cpu_wrapper_chunked_fast)
    drift_optimizer.LBFGSB_OPTIONS = dict(drift_optimizer.LBFGSB_OPTIONS,
                                          ftol=float(settings.optimizer_ftol))
    drift_optimizer.LBFGSB_OPTIONS_COARSE = None if settings.optimizer_ftol_coarse is None else dict(
        drift_optimizer.LBFGSB_OPTIONS, ftol=float(settings.optimizer_ftol_coarse))

    # the drift must cover every frame of the table it will be applied to,
    # whatever subset it is estimated from
    all_frames = np.asarray(locs["frame"], dtype=np.int64)
    n_frames = int(all_frames.max()) + 1 if len(all_frames) else 0

    backend = settings.backend or best_backend()
    if settings.quality_control:
        backend += "_qc"
    drift_optimizer.LAST_QUALITY_CONTROL = None
    drift_optimizer.FLAG_SEGMENTS = lambda q_obs, q_null: flag_flawed_segments_by_lift(
        q_obs, q_null, settings.min_lift)

    source = locs[_selection(locs, select)]
    if settings.group:
        source = _grouped(source, settings, pixelsize_nm)
    x, y, z = _positions_nm(source, pixelsize_nm, settings.use_z)
    dataset = np.column_stack([x, y, z,
                               np.asarray(source["frame"], dtype=np.float64)])
    if len(dataset) < 2:
        raise ValueError("need at least two localizations to estimate drift")

    # COMET's `display` is both its progress log and a blocking `plt.show()`;
    # we want the log, and to draw the drift ourselves when asked (Drift.plot)
    with _without_blocking_plots(display):
        drift = comet_run_kd(
            dataset,
            segmentation_mode=settings.segmentation_mode,
            segmentation_var=settings.segmentation_var,
            max_drift_nm=settings.max_drift_nm,
            target_sigma_nm=settings.target_sigma_nm,
            initial_sigma_nm=settings.initial_sigma_nm,
            boxcar_width=settings.boxcar_width,
            interpolation_method=settings.interpolation,
            max_locs_per_segment=settings.max_locs_per_segment,
            mode=backend,
            display=display,
            # an upfront pair-count estimate: the neighbour search runs over the
            # whole dataset at once, and is where too large a max_drift_nm turns
            # into minutes and gigabytes.  Better warned before than after.
            pair_indices_safety_check=True,
            # the estimate is made from a subset, but must cover every frame of
            # the table it will be applied to
            min_max_frames=(0, n_frames - 1),
        )
    quality = drift_optimizer.LAST_QUALITY_CONTROL
    return Drift(drift[:, :3], settings, n_used=len(dataset),
                 flagged_windows=None if quality is None else quality["flagged"])


def _spline_basis(x, n_coefficients: int, lo: float, hi: float, degree: int = 3):
    """Cubic B-spline design matrix: row per point, column per coefficient."""
    from scipy.interpolate import BSpline

    interior = np.linspace(lo, hi, n_coefficients - degree + 1)
    knots = np.concatenate([np.full(degree, lo), interior, np.full(degree, hi)])
    return BSpline.design_matrix(np.clip(np.asarray(x, dtype=float), lo, hi),
                                 knots, degree)


def _estimate_spline(locs: Localizations, settings: DriftSettings, select,
                     pixelsize_nm: Optional[float], display: bool) -> "Drift":
    """Fit drift as a spline in time, with no segmentation at all.

    The optimizer's variables are the spline's coefficients; the cost and its
    gradient are COMET's, over the same neighbour pairs.  Each localization
    carries its own frame, so `mu` handed to the kernel is the drift *per
    frame*, and the chain rule to the coefficients is one matrix product.
    """
    from comet.core.cpu_wrapper import cpu_wrapper_chunked_approx
    from comet.core.pair_indices import pair_indices_kdtree
    from scipy.optimize import minimize

    all_frames = np.asarray(locs["frame"], dtype=np.int64)
    n_frames = int(all_frames.max()) + 1 if len(all_frames) else 0

    source = locs[_selection(locs, select)]
    if settings.group:
        source = _grouped(source, settings, pixelsize_nm)
    x, y, z = _positions_nm(source, pixelsize_nm, settings.use_z)
    coords = np.ascontiguousarray(np.column_stack([x, y, z]), dtype=np.float32)
    frames = np.ascontiguousarray(source["frame"], dtype=np.int32)

    idx_i, idx_j, ok = pair_indices_kdtree(coords, settings.max_drift_nm)
    if not ok:
        raise MemoryError("the neighbour search ran out of memory; reduce "
                          "max_drift_nm or use group=True")

    n_coefficients = max(4, int(round(n_frames / settings.spline_knot_frames)) + 3)
    basis = _spline_basis(np.arange(n_frames), n_coefficients, 0.0, n_frames - 1.0)
    basis_t = basis.T.tocsr()
    _log(display, f"spline drift: {n_coefficients} coefficients over {n_frames} "
                  f"frames, {len(idx_i):,} pairs")

    penalty = float(settings.spline_penalty)
    if penalty:
        # second differences of the coefficients: what "not oscillating" means
        d2 = np.diff(np.eye(n_coefficients), 2, axis=0)
        curvature = d2.T @ d2

    def objective(flat, sigma, factor):
        c = flat.reshape(n_coefficients, 3)
        mu = np.ascontiguousarray(basis @ c)
        value, gradient = cpu_wrapper_chunked_approx(
            mu, coords, frames, idx_i, idx_j, sigma, factor)
        gradient = basis_t @ gradient.reshape(n_frames, 3)
        if penalty:
            value += penalty * float((c * (curvature @ c)).sum())
            gradient = gradient + 2.0 * penalty * (curvature @ c)
        return value, gradient.ravel()

    sigma = (settings.initial_sigma_nm if settings.initial_sigma_nm is not None
             else settings.max_drift_nm / 3.0)
    target = settings.target_sigma_nm
    bound = settings.max_drift_nm * 2.0
    coefficients = np.zeros(n_coefficients * 3)
    while True:
        result = minimize(objective, coefficients, args=(sigma, 1.0), jac=True,
                          method="L-BFGS-B", bounds=[(-bound, bound)] * coefficients.size,
                          options={"ftol": settings.optimizer_ftol, "gtol": 1e-5,
                                   "maxls": 40})
        coefficients = result.x
        _log(display, f"  sigma {sigma:6.1f} nm: {result.nfev} evaluations, "
                      f"cost {result.fun:.6g}")
        if sigma <= target:
            break
        sigma = max(sigma / 1.5, target)

    drift = np.asarray(basis @ coefficients.reshape(n_coefficients, 3))
    return Drift(drift, settings, n_used=len(coords))


def _log(enabled: bool, message: str) -> None:
    if enabled:
        print(message)


def _two_stage(locs: Localizations, settings: DriftSettings, select,
               pixelsize_nm: Optional[float], display: bool) -> "Drift":
    """A grouped pass to get close, then an ungrouped one over a small radius."""
    coarse = estimate_drift(locs, replace(settings, two_stage=False, group=True),
                            select, pixelsize_nm, display)
    radius = settings.two_stage_radius_nm
    fine = estimate_drift(
        coarse.apply(locs, pixelsize_nm),
        replace(settings, two_stage=False, group=False, max_drift_nm=radius,
                initial_sigma_nm=radius / 3.0, target_sigma_nm=radius / 5.0),
        select, pixelsize_nm, display)
    # the two are drift of the same sample measured in sequence, so they add
    return Drift(coarse.drift + fine.drift, settings, n_used=fine.n_used,
                 flagged_windows=fine.flagged_windows)


def correct_drift(locs: Localizations, settings: Optional[DriftSettings] = None,
                  select=None, pixelsize_nm: Optional[float] = None,
                  display: bool = False):
    """Estimate from the selection, correct everything.  Returns (locs, drift)."""
    drift = estimate_drift(locs, settings, select, pixelsize_nm, display)
    return drift.apply(locs, pixelsize_nm), drift


# ----------------------------------------------------------------------- saving
def drift_corrected_path(path) -> Path:
    """``foo.hdf5`` -> ``foo_driftc.hdf5``."""
    path = Path(path)
    return path.with_name(path.stem + "_driftc" + (path.suffix or ".hdf5"))


def save_drift_corrected(path, locs: Localizations, drift: Drift) -> Path:
    """Write the corrected table, with the drift curve alongside it.

    The file is an ordinary smapfit localization file -- the viewer and every
    other reader take it as it is -- plus a ``/drift`` group, so the correction
    that was applied can be inspected or undone later.
    """
    import h5py

    from .io.hdf5 import save_localizations

    path = save_localizations(path, locs)
    with h5py.File(path, "a") as f:
        group = f.create_group("drift")
        for axis, name in enumerate(("x_nm", "y_nm", "z_nm")):
            group.create_dataset(name, data=drift.drift[:, axis].astype(np.float32))
        group.create_dataset("frame", data=drift.frames.astype(np.int64))
    return path


def load_drift(path) -> Drift:
    """Read the drift curve stored next to a corrected table."""
    import h5py

    with h5py.File(path, "r") as f:
        group = f["drift"]
        drift = np.column_stack([group[name][()]
                                 for name in ("x_nm", "y_nm", "z_nm")])
    return Drift(drift)


# ---------------------------------------------------------------------- helpers
def _grouped(locs: Localizations, settings: "DriftSettings",
             pixelsize_nm: Optional[float]) -> Localizations:
    """Collapse each emitter's blink into one localization before estimating."""
    from .group import GroupSettings, group

    dx = settings.group_dx_nm
    if "x_nm" not in locs:  # a table in pixels; the linking radius is in nm
        dx /= _pixelsize_nm(locs, pixelsize_nm)
    grouped, _ = group(locs, GroupSettings(dx=dx, dt=settings.group_dt))
    return grouped


@contextmanager
def _without_blocking_plots(active: bool):
    """Let COMET log its progress without stopping in its own `plt.show()`.

    Its diagnostic figures are closed again; a drift curve is one line here
    (`Drift.plot`) and belongs where the caller decides to draw it.
    """
    if not active:
        yield
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        yield
        return
    existing = set(plt.get_fignums())
    show = plt.show
    plt.show = lambda *args, **kwargs: None
    try:
        yield
    finally:
        plt.show = show
        for number in set(plt.get_fignums()) - existing:
            plt.close(number)


def _selection(locs: Localizations, select) -> np.ndarray:
    """A LocFilter, a boolean mask or an index array -> indices."""
    n = len(locs)
    if select is None:
        return np.arange(n)
    indices = getattr(select, "indices", None)  # a LocFilter
    if indices is not None:
        return np.asarray(indices)
    select = np.asarray(select)
    if select.dtype == bool:
        if select.shape != (n,):
            raise ValueError(f"mask has {select.shape[0]} entries, table has {n}")
        return np.flatnonzero(select)
    return select.astype(np.int64)


def _pixelsize_nm(locs: Localizations, pixelsize_nm: Optional[float]) -> float:
    """The pixel size, given or from the table's own provenance."""
    if pixelsize_nm is not None:
        return float(pixelsize_nm)
    value = locs.metadata.get("pixelsize_nm")
    if value is None:
        camera = locs.metadata.get("camera") or {}
        micron = camera.get("pixelsize_um") if isinstance(camera, dict) else None
        value = None if micron is None else float(micron) * 1000.0
    if value is None:
        raise ValueError("this table is in pixels and its pixel size is not "
                         "recorded; pass pixelsize_nm=...")
    return float(value)


def _positions_nm(locs: Localizations, pixelsize_nm: Optional[float],
                  use_z: Optional[bool]):
    """x, y, z of the whole table in nm; z is zero for a 2D table."""
    if "x_nm" in locs and "y_nm" in locs:
        x = np.asarray(locs["x_nm"], dtype=np.float64)
        y = np.asarray(locs["y_nm"], dtype=np.float64)
    elif "x_pix" in locs and "y_pix" in locs:
        scale = _pixelsize_nm(locs, pixelsize_nm)
        x = np.asarray(locs["x_pix"], dtype=np.float64) * scale
        y = np.asarray(locs["y_pix"], dtype=np.float64) * scale
    else:
        raise ValueError("no position columns (x_nm/y_nm or x_pix/y_pix)")

    has_z = "z_nm" in locs
    if use_z is None:
        use_z = has_z
    elif use_z and not has_z:
        raise ValueError("use_z=True but the table has no z_nm column")
    z = np.asarray(locs["z_nm"], dtype=np.float64) if use_z else np.zeros_like(x)
    return x, y, z
