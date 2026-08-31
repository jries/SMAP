import argparse
from comet.core.backends import best_backend
from comet.core.drift_optimizer import comet_run_kd
from comet.core.io_utils import (
    load_thunderstorm_csv,
    load_normal_molecule_set,
    correct_and_save_thunderstorm_csv,
    save_dataset_as_thunderstorm_csv,
)


def main():
    parser = argparse.ArgumentParser(description="COMET drift correction CLI")

    parser.add_argument("--input", "-i", required=True, help="Path to input file (.csv or .h5)")
    parser.add_argument("--output", "-o",  required=True, help="Path to save corrected localizations")
    parser.add_argument("--segmentation_mode", "-sm", type=int, default=2, help="0: num windows, 1: locs/window, "
                                                                         "2: frames/window")
    parser.add_argument("--segmentation_var", "-sv", type=int, required=True, help="Segmentation variable (depends on mode)")
    parser.add_argument("--initial_sigma_nm", "-isig", type=float, default=None,
                        help="Initial sigma (nm). Defaults to max_drift_nm/3 when omitted.")
    parser.add_argument("--target_sigma_nm", "-tsig", type=float, default=30.0)
    parser.add_argument("--max_drift_nm", "-d", type=float, default=300.0,
                        help="Maximum expected drift (nm); also the neighbour-search radius")
    parser.add_argument("--boxcar_width", type=int, default=1)
    parser.add_argument("--interpolation", choices=["cubic", "catmull-rom"], default="cubic")
    parser.add_argument("--format", choices=["csv", "h5"], help="Output file format", default="csv")
    parser.add_argument("--pixelsize_nm", type=float, default=160.0,
                        help="Camera pixel size in nm, recorded in --format h5 output "
                             "(molecule sets store positions in pixel units)")
    parser.add_argument("--pixelsize_z_nm", type=float, default=None,
                        help="Axial pixel size in nm for --format h5 output (default: same as --pixelsize_nm)")
    parser.add_argument("--display", action="store_true", help="Display intermediate results during processing")
    parser.add_argument("--max_locs_per_segment", "-mlpseg", type=int, default=None)
    parser.add_argument("--mode", choices=["cuda", "cuda_qc", "cpu", "torch", "torch_qc"], default=None,
                        help="Compute backend (default: fastest available on this machine)")
    args = parser.parse_args()

    if args.mode is None:
        args.mode = best_backend()
        if args.display:
            print(f"Using '{args.mode}' backend.")

    input_is_csv = args.input.lower().endswith(".csv")
    input_is_h5 = args.input.lower().endswith(".h5")
    if input_is_csv:
        dataset = load_thunderstorm_csv(args.input)
    elif input_is_h5:
        dataset = load_normal_molecule_set(args.input)
    else:
        raise ValueError("Unsupported input format; expected .csv or .h5.")

    drift, corrected = comet_run_kd(
        dataset=dataset,
        segmentation_mode=args.segmentation_mode,
        segmentation_var=args.segmentation_var,
        initial_sigma_nm=args.initial_sigma_nm,
        target_sigma_nm=args.target_sigma_nm,
        max_drift_nm=args.max_drift_nm,
        max_locs_per_segment=args.max_locs_per_segment,
        boxcar_width=args.boxcar_width,
        return_corrected_locs=True,
        interpolation_method=args.interpolation,
        mode=args.mode,
        display=args.display,
        save_corrected_locs=args.format == "h5",
        # without this the h5 path ignores --output and opens a save dialog
        save_filepath=args.output if args.format == "h5" else None,
        pixelsize_nm=args.pixelsize_nm,
        pixelsize_z_nm=args.pixelsize_z_nm,
    )

    if args.format == "csv":
        if input_is_csv:
            # Re-reads the source CSV so every original column is carried through,
            # not just x/y/z/frame.
            correct_and_save_thunderstorm_csv(drift, args.input, args.output)
        else:
            # An h5 input has no CSV to re-read; write the corrected array directly.
            save_dataset_as_thunderstorm_csv(corrected, args.output)


if __name__ == "__main__":
    main()
