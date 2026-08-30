# Design notes

Decisions that are not obvious from the code, and open questions.  Everything
here was decided deliberately; if something looks wrong, read this first.

## Scope

A stripped-down version of SMAP's `fit_fastsimple` workflow: TIFF -> photons ->
filter -> peak finding -> ROI cutting -> MLE fit -> HDF5.  Single channel,
Poisson noise only.  Fit models: `gauss_free`, `gauss_xy`, `cspline`.

Deliberately **not** ported: sCMOS per-pixel variance, the Anscombe transform
and probability cutoff, background estimation, minimum-distance filtering,
mean-PSF / MIP-PSF detection filters, ROI masks, multi-channel / global fitting,
Zernike and astigmatic-Gaussian z models, multiple z start values.

## Departures from SMAP, and why

* **No module chain.**  SMAP's WorkflowModules exist to serve a GUI (parameter
  synchronisation, per-frame data packets, global-variable accumulators).  Here
  the stages are plain functions over blocks of frames, and the only state that
  survives between blocks -- the ROI buffer -- lives in `LocalizationEngine`.
* **Chunked, not per-frame.**  Every stage takes `(n, y, x)`.  This is what makes
  the vectorised and threaded implementations possible.
* **No x/y clamp in the fitter.**  SMAP's single-channel kernels clamp x and y to
  the central half of the ROI; the multi-channel fitter dropped that and we
  follow it.  A runaway fit is now visible instead of silently parked on a
  boundary -- reject them with `FitSettings.max_fit_distance` if wanted.
* **Border candidates are dropped, not clamped.**  SMAP shifts a ROI inward and
  fits it anyway, which biases those localizations because the emitter is no
  longer near the ROI centre.
* **No background offset before filtering.**  In SMAP that value (the minimum of
  the image's first column) only compensated the zero padding of the
  convolution at the borders; we pad by edge replication instead.  Verified: the
  interior is unaffected to 2e-5, only a ~3 px rim differs.
* **EM excess noise lives outside the fitter.**  Photons are divided by 2 before
  fitting and photons/background scaled back afterwards, so the fitter assumes
  pure Poisson statistics and knows nothing about cameras.
* **Pixels are primary.**  The fit measures pixels; nm is derived at the very end
  via `to_nm()` / `FitSettings.output_unit`, because the pixel size is a separate
  calibration that may be corrected later.  `z_nm` is the exception: it comes
  from the calibration's `dz` and has no pixel equivalent.
* **No mirroring in the fitter.**  A calibration built from mirrored bead images
  (`parameters.emmirror` in the `_3Dcal.mat`) is handled by flipping the ROI in x
  and flipping the fitted x back; the fitter itself is orientation-free.

## Conventions

Fixed once, in `io/calibration.py`, so nothing downstream deals with MATLAB's
layout:

* images `(y, x)` C-contiguous; spline coefficients `(64, nz, ny, nx)`, x fastest
* fit parameters `(x, y, photons, background, [z | sigma | sigma_x, sigma_y])`
* MATLAB's first spatial axis is the image **row**, which is why `P(:,1)` is y
  there; both the spatial axes and the 64 monomial indices are swapped on load

## Verification that is not circular

* the loaded spline coefficients reproduce the bead stack stored beside them in
  the same file (correlation 0.9999, vs 0.82/0.82/0.83 for transposed,
  x-mirrored and z-flipped controls) -- `scripts/check_calibration.py`
* the Gaussian models are checked against independently simulated data: unbiased
  to <0.01 px, scatter matches the returned CRLB
* the mirror flip was confirmed on real data: the flipped orientation has the
  better log-likelihood for 65.8% of ROIs (median dlogL +2.82)
* threaded output is bit-identical to serial, at every stage

## Performance

M1 Pro, 46,005 frames of 200x200, spline fit, ROI 13: **35.9 s** (1283 frames/s,
~850k localizations).  Detection 12.9 s, fitting 22.2 s, I/O 0.8 s.

Findings worth keeping:

* the C++ single-pass maximum search beats `scipy.ndimage.maximum_filter` 5.2x
  (short-circuit evaluation rejects most pixels after one comparison)
* the 2-D DoG kernel has **rank 2**, so it is not separable -- subtracting the
  1-D kernels is wrong.  Subtracting the 2-D kernels and doing one convolution
  is exact and, below radius 4, faster than the fused separable passes
  (`_SEPARABLE_FROM_RADIUS`)
* SMAP's DoG window is `max(ceil(6*sigma-1), 3)`: radius 3 for sigma 1.2
* threading gives ~4x on filtering and maxima, ~6.5x on fitting
* the per-frame dynamic-cutoff loop is the last serial part (~2.7 s, 80% of it
  numpy dispatch overhead).  It could be vectorised with one lexsort per chunk.

## Open questions

* Fitted x sits ~0.24 px from the peak-finder position, and the sign flips with
  the mirror (-0.22 unmirrored, +0.26 mirrored) while y stays at -0.09.  A
  round trip (render the model, fit it back) is exact, so this is a property of
  the data/calibration, not of the code.  Worth comparing against SMAP's own
  output on the same dataset.
* Camera offset: taken from the image metadata, but an iXon does not report one,
  so the settings file is used with a warning.  Should that be a hard error?
* A configuration layer that merges file metadata with user settings once, up
  front, so downstream code never deals with partial metadata.

## Environment

`python3` on this machine is x86_64 conda under Rosetta.  Use the native arm64
venv:

    .venv/bin/python -m pytest tests/
    .venv/bin/python setup.py build_ext --build-lib src

The extension builds universal2.  `setup.py` contains a workaround for a
Command Line Tools install with incomplete libc++ headers; it is gated on a
test-compile and stays silent on a healthy toolchain.
