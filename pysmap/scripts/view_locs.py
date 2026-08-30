#!/usr/bin/env python3
"""Open the interactive viewer on a saved localization file.

    view_locs.py OUT.h5 [--mode precision|gauss|hist] [--color z_nm] [--lut hot]

Scroll or pinch to zoom about the cursor, +/- to zoom about the centre, drag to
pan, r to reset.  The sliders filter; a slider at either end means no bound.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from smapfit.io.hdf5 import load_localizations
from smapfit.render import DisplaySettings, RenderSettings
from smapfit.viewer import show


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file", help="localizations saved by fit_dataset.py")
    p.add_argument("--mode", default="precision",
                   choices=["hist", "gauss", "precision"])
    p.add_argument("--sigma", type=float, default=10.0,
                   help="rendering sigma for --mode gauss, in the file's units")
    p.add_argument("--color", default=None,
                   help="colour by this field instead of by density, e.g. z_nm")
    p.add_argument("--range", type=float, nargs=2, default=None,
                   metavar=("MIN", "MAX"), help="colour range for --color")
    p.add_argument("--lut", default=None, help="colour map (default: hot, or "
                                               "turbo when colouring by a field)")
    p.add_argument("--additive", action="store_true",
                   help="SMAP's colour composite: overlapping colours add, so "
                        "red over cyan saturates to white")
    p.add_argument("--contrast", type=float, default=None,
                   help="dynamic contrast p: saturates 10^-p of the pixels "
                        "(default 3; lower is brighter)")
    args = p.parse_args()

    locs = load_localizations(args.file)
    print(locs)
    print(f"units: {locs.metadata.get('units', 'unknown')}")

    settings = RenderSettings(mode=args.mode, sigma=args.sigma,
                              color_field=args.color, color_range=args.range)
    lut = args.lut or ("turbo" if args.color else "hot")
    display = DisplaySettings(lut=lut, color_mode="sum" if args.additive else "hue")
    if args.contrast is not None:
        display.contrast = args.contrast
    show(locs, settings, display)


if __name__ == "__main__":
    main()
