#!/usr/bin/env python3
"""Drift-correct a saved localization file with COMET.

    drift_correct.py OUT.hdf5 [--filter FIELD LO HI]... [--frames-per-window 500]

The drift is estimated from the localizations that pass the filters -- give the
same limits you would use in the viewer -- and then subtracted from *all* of
them.  The result is written next to the input as ``..._driftc.hdf5``, with the
drift curve stored in the file.

Needs COMET:  pip install -e externaltools/Comet/Python_interface
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from smappy.drift import DriftSettings, correct_drift, drift_corrected_path, save_drift_corrected
from smappy.rcc import RCCSettings, estimate_drift_rcc
from smappy.filter import LocFilter
from smappy.io.hdf5 import load_localizations


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file", help="localizations saved by fit_dataset.py")
    p.add_argument("--out", default=None, help="output file (default: ..._driftc.hdf5)")
    p.add_argument("--filter", nargs=3, action="append", default=[],
                   metavar=("FIELD", "MIN", "MAX"),
                   help="estimate only from localizations in this range; "
                        "'-' is unbounded.  Repeatable, e.g. "
                        "--filter loc_precision_nm - 20 --filter logl_rel -2 -")

    windows = p.add_mutually_exclusive_group()
    windows.add_argument("--frames-per-window", type=int, default=500,
                         help="time window size in frames (default 500)")
    windows.add_argument("--locs-per-window", type=int, default=None,
                         help="time window size in localizations instead")
    windows.add_argument("--windows", type=int, default=None,
                         help="fixed number of time windows instead")

    p.add_argument("--max-drift", type=float, default=300.0,
                   help="largest expected drift in nm, also the neighbour "
                        "search radius (default 300)")
    p.add_argument("--target-sigma", type=float, default=30.0,
                   help="length scale the refinement stops at, in nm (default 30)")
    p.add_argument("--smooth", type=int, default=1,
                   help="boxcar width in time windows (default 1, no smoothing)")
    p.add_argument("--max-locs-per-window", type=int, default=None,
                   help="cap per window, to bound memory and runtime")
    p.add_argument("--no-z", action="store_true",
                   help="estimate lateral drift only, even for a 3D dataset")
    p.add_argument("--backend", default=None, choices=["cuda", "torch", "cpu"],
                   help="default: the fastest available")
    p.add_argument("--pixelsize", type=float, default=None,
                   help="nm per pixel, if the file is in pixels and does not "
                        "record it")
    p.add_argument("--ftol", type=float, default=DriftSettings.optimizer_ftol,
                   help="L-BFGS-B relative tolerance (default %(default)g; "
                        "COMET's own 2e-13 is ~2.7x slower for a sub-nm change)")
    p.add_argument("--rcc", action="store_true",
                   help="estimate by redundant cross-correlation of rendered "
                        "time windows instead of COMET; an independent method, "
                        "useful as a second opinion on a drift curve")
    p.add_argument("--rcc-windows", type=int, default=RCCSettings.n_timepoints,
                   help="time windows for --rcc (default %(default)d)")
    p.add_argument("--per-window", action="store_true",
                   help="fit one free drift vector per time window instead of "
                        "the default spline in time")
    p.add_argument("--ungrouped", action="store_true",
                   help="estimate from every localization instead of one per "
                        "blink (slower, and measurably noisier)")
    p.add_argument("--knot-frames", type=int, default=DriftSettings.spline_knot_frames,
                   help="spline coefficient spacing in frames (default %(default)d)")
    p.add_argument("--two-stage", action="store_true",
                   help="grouped pass, then an ungrouped one over a small "
                        "radius: ~9x faster for ~99%% of the improvement")
    p.add_argument("--qc", action="store_true",
                   help="discard time windows whose estimate does not beat the "
                        "no-correction overlap, and interpolate across them")
    p.add_argument("--plot", action="store_true", help="show the drift curve")
    a = p.parse_args()

    locs = load_localizations(a.file)
    print(locs)

    ranges = {field: (_bound(lo), _bound(hi)) for field, lo, hi in a.filter}
    keep = LocFilter(locs, **ranges) if ranges else None
    if keep is not None:
        print(f"estimating from {len(keep.indices)} of {len(locs)} localizations "
              + ", ".join(f"{f} in [{lo}, {hi}]" for f, (lo, hi) in ranges.items()))

    if a.windows is not None:
        mode, var = 0, a.windows
    elif a.locs_per_window is not None:
        mode, var = 1, a.locs_per_window
    else:
        mode, var = 2, a.frames_per_window

    settings = DriftSettings(
        segmentation_mode=mode, segmentation_var=var,
        max_drift_nm=a.max_drift, target_sigma_nm=a.target_sigma,
        boxcar_width=a.smooth, max_locs_per_segment=a.max_locs_per_window,
        backend=a.backend, use_z=False if a.no_z else None,
        optimizer_ftol=a.ftol, group=not a.ungrouped, quality_control=a.qc,
        two_stage=a.two_stage, spline=not a.per_window,
        spline_knot_frames=a.knot_frames)

    if a.rcc:
        drift = estimate_drift_rcc(
            locs, RCCSettings(n_timepoints=a.rcc_windows,
                              use_z=False if a.no_z else None),
            select=keep, pixelsize_nm=a.pixelsize, display=True)
        corrected = drift.apply(locs, a.pixelsize)
    else:
        corrected, drift = correct_drift(locs, settings, select=keep,
                                         pixelsize_nm=a.pixelsize, display=True)
    print(drift)

    out = save_drift_corrected(a.out or drift_corrected_path(a.file), corrected, drift)
    print(f"wrote {out}")

    if a.plot:
        import matplotlib.pyplot as plt
        drift.plot()
        plt.show()


def _bound(value):
    return None if value in ("-", "", "none") else float(value)


if __name__ == "__main__":
    main()
