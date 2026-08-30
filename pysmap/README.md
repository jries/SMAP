# smapfit

Python implementation of the SMAP single-molecule fitting pipeline: camera
conversion, filtering, peak finding, ROI cutting, maximum-likelihood fitting
with a Gaussian or experimental (cubic-spline) PSF, and streaming output to
HDF5.  Reads SMAP `_3Dcal.mat` calibration files and Micro-Manager TIFF stacks.

See [NOTES.md](NOTES.md) for the design decisions and open questions.

## Install

    /usr/bin/python3 -m venv .venv                 # native arm64 on Apple silicon
    .venv/bin/python -m pip install numpy scipy tifffile h5py pybind11 pyyaml pytest
    .venv/bin/python setup.py build_ext --build-lib src

## Use

    from smapfit.io.tiff import open_stack, camera_metadata
    from smapfit.io.cameras_mat import CameraPresets
    from smapfit.io.calibration import load_spline_calibration
    from smapfit.detect import DoGFilter, DynamicCutoff, PeakFinder
    from smapfit.psf import SplinePSF, GaussianPSF
    from smapfit.pipeline import FitSettings, fit_stack

    source = open_stack("...MMStack_Default.ome.tif")
    camera = camera_metadata(source, CameraPresets.load("RiesLab_cameras.mat"),
                             {"pixelsize_um": 0.127})
    model = SplinePSF(load_spline_calibration("..._3dcal.mat"))
    finder = PeakFinder(DoGFilter(1.2), DynamicCutoff(1.7))

    locs, engine = fit_stack(source.frames(chunk=200), camera, finder, model,
                             FitSettings(roisize=13, output_unit="nm"))

Or from the command line:

    .venv/bin/python scripts/fit_dataset.py DATA CAMERAS.mat CONFIG.yaml OUT.h5 \
        --cal CAL_3dcal.mat --units nm

For online analysis, drive `LocalizationEngine` directly: `push(frames)` returns
localizations once enough ROIs have accumulated, `flush()` forces a partial
block.  Nothing asks how many frames there will be.

## Scripts

| script | what it checks |
|---|---|
| `check_calibration.py` | spline coefficients against the bead stack in the same file |
| `check_stack.py` | what the image metadata provides, and what is missing |
| `check_detection.py` | filtering, peak finding and ROI cutting on real frames |
| `check_fit.py` | spline fits on real data, with and without the mirror flip |
| `fit_dataset.py` | the whole pipeline, to HDF5 |

## Tests

    SMAPFIT_TEST_CAL=/path/to/_3dcal.mat PYTHONPATH=src .venv/bin/python -m pytest tests/
