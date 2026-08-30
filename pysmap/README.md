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

To render an image from a localization table:

    from smapfit.filter import LocFilter
    from smapfit.render import (FieldOfView, RenderSettings, DisplaySettings,
                                render_locs)

    keep = LocFilter(locs, loc_precision_nm=(None, 20), logl_rel=(-2, 0))
    fov = FieldOfView.around(locs["x_nm"], locs["y_nm"], pixelsize=10.0)
    image = render_locs(locs, fov, RenderSettings(mode="precision"), select=keep)
    rgb = DisplaySettings(lut="hot", gamma=0.7).apply(image)

`mode` is `"hist"`, `"gauss"` (one sigma for all) or `"precision"` (sigma from
the localization precision, the default in SMAP).  Set `color_field` to colour
by z or any other column instead of by density.  Rendering and display are
separate on purpose: contrast, gamma and the colour map change without
re-rendering.

To merge localizations of the same emitter across consecutive frames:

    from smapfit.group import group, GroupSettings

    grouped, group_index = group(locs, GroupSettings(dx=50.0, dt=1))

`grouped` carries the same columns, combined by SMAP's per-column rules
(positions weighted by precision, z by its own error, photons summed and their
errors added in quadrature, precisions added in inverse quadrature), plus
`n_in_group`.

To look at the result:

    from smapfit.viewer import show
    show(locs)                       # or: scripts/view_locs.py OUT.h5

Scroll or pinch to zoom about the cursor, `+`/`-` to zoom about the centre, drag
to pan, `r` to reset.  Type a minimum and maximum to filter on localization
precision, z, PSF size, relative log-likelihood and frame; an empty box means
"no bound", and the data range is shown beside each row.
The "grouped" box switches to the grouped table, which is built on first use and
keeps its own filter; "additive" switches field colouring to SMAP's composite,
where overlapping colours add (red over cyan saturates to white).

"colour by" selects a plain intensity image or one coded by z, frame,
localization precision or photons, with the range typed into the "colour" row
(empty ends fall back to the data's own).  The range is always explicit, so the
same z means the same colour at every zoom and after every block of a live fit;
the LUT follows the choice -- `hot` for intensity, `turbo` for a coded field.
Needs matplotlib (`pip install matplotlib`).

Or from the command line:

    .venv/bin/python scripts/fit_dataset.py DATA CAMERAS.mat CONFIG.yaml OUT.h5 \
        --cal CAL_3dcal.mat --units nm

## Online: fit while the microscope writes

    .venv/bin/python scripts/live_fit.py DATA CAMERAS.mat CONFIG.yaml OUT.h5 \
        --cal CAL_3dcal.mat --update 3 --timeout 30

`DATA` is the growing Micro-Manager TIFF, or the directory it is being written
into; it does not have to exist yet.  The window opens as soon as the first
frames appear and takes in new localizations every `--update` seconds; the fit
ends `--timeout` seconds after the last frame is written, which is how an
acquisition stops.  `OUT.h5` is written throughout and is the result.

Everything the offline viewer offers works while this runs -- zoom, pan, the
filter boxes, contrast, grouping -- and **an update changes none of them**: new
localizations appear inside the view being looked at, under the bounds already
typed.  The frame comes from the camera field of view, so the image does not
rescale as data arrives.  Grouping cannot be extended, so the grouped table is
marked stale and rebuilt when it is next asked for.

From Python:

    from smapfit.live import LiveSettings, live_view

    live_view(directory, camera, finder, model, FitSettings(output_unit="nm"),
              output="OUT.h5", live=LiveSettings(update_seconds=3.0))

`LiveFit` is the same thing without a window: it runs the pipeline in a thread
and queues finished blocks, for a different front end or a headless run.

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
| `view_locs.py` | opens the viewer on a saved localization file |

## Tests

    SMAPFIT_TEST_CAL=/path/to/_3dcal.mat PYTHONPATH=src .venv/bin/python -m pytest tests/
