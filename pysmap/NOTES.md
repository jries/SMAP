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

## Rendering

The renderer is deliberately two stages: accumulation into float planes, and a
display step (normalise, gamma, LUT) that turns them into RGB.  Changing the
contrast or the colour map therefore never re-renders, which is what a viewer
needs.  SMAP splits the same way (`renderSMAP` / `drawerSMAP`).

* **One kernel for all three modes.**  The kernel is the *pixel integral* of the
  Gaussian, `erf(right) - erf(left)`, not a point sample of it.  A rendered
  pixel is an integral over its area, and point sampling breaks down exactly
  where sigma approaches one pixel -- which is where the renderer spends most of
  its time.  It also makes the histogram the sigma -> 0 limit of the Gaussian
  rather than a separate model, and both are tested against each other.
* **No template.**  SMAP looks up a 601x601 Gaussian template with
  nearest-neighbour sampling.  The kernel is separable, so a (2d+1)^2 ROI costs
  2(2d+2) `erf` calls plus the outer product, which dominates either way: the
  template buys nothing and costs quantisation error.
* **Normalised by its own sum, not by an erf truncation correction.**  Every
  localization then contributes exactly N, whatever its sigma or subpixel
  position, so histogram and Gaussian images carry the same total intensity.
  The kernel is built over the full ROI even where the image clips it, so a
  localization at the border loses the outside part instead of having it
  redistributed inwards -- and the result does not depend on thread boundaries.
* **Field colouring keeps four planes**, the three colour-weighted ones plus the
  plain weight, so both composites come from one render and the switch between
  them is a re-display.  `sum` is SMAP's: display the colour planes directly, so
  a red and a cyan localization in one saturated pixel **add to white**, as
  additive colour should.  `hue` divides by the weight first, so that pixel
  stays a mid grey at full brightness and the hue never drifts towards white
  with density.  Below saturation the two are algebraically identical -- `hue`
  divides by the weight and multiplies it back in -- so the difference is
  confined to the brightest pixels: on the clathrin data at p = 2.5, 4-8% of
  pixels differ by more than 0.02, and those are the dense cores and the
  fiducials.  Which is right depends on what the image is for: `sum` reads as
  brightness-is-density with colour bleaching out where it is densest, `hue`
  keeps z legible in exactly those cores.  `hue` is the default; the viewer's
  "additive" box switches.
* **Contrast is SMAP's dynamic contrast**, parameterised the way SMAP's
  `imax_min` is: one number p saturates 10^-p of the pixels.  The value
  distribution of a superresolution image is extremely heavy-tailed -- on real
  data a fiducial bead outweighs the structure by three orders of magnitude --
  so an absolute maximum does not transfer between datasets and the useful
  range of p spans decades, which is why the slider moves the exponent.  p = 3
  is the default; SMAP's own scripts use 3.5, which is too dark for the
  clathrin dataset.
* **The pixel grid is stated once.**  A position lands in pixel
  `floor((x - x0) / pixelsize)`.  SMAP rounds and shifts the range by half a
  pixel first, which is the same grid expressed twice.
* Not ported: the transparency/occlusion modes (`gaussrenderT`,
  `gaussrenderTx`), layers and compositing, the DL/tiff modes, grouped
  localizations, `normalizeFoV`, `remout`.  The sigma policy (`gaussfac`,
  `mingaussnm`, `mingausspix`, the cap at 10x the median) is kept, in
  `SigmaSettings`, because it is what makes precision-weighted images readable.

**Filtering** (`LocFilter`) caches one boolean array per field and recomputes
only the field that changed, which is what SMAP does and the reason interactive
filtering is possible at all; the displayed set is an AND over a handful of
arrays.  `indices` hands the renderer a compacted selection so filtering indexes
three or four columns instead of copying the table.  A range excludes NaN,
because a comparison against NaN is false -- a localization with no z should not
appear in a z slab.  Not ported: ROI/polygon masks and the grouped-localization
filters; `set_mask` takes an arbitrary boolean array under a name, which is
where those, or the viewer's "inside the field of view", would plug in.

Changing the LUT is free -- it is applied to the finished image -- *except* when
colouring by a field, where the colour is baked into the accumulation and the
image has to be rendered again.  That is why `render_locs` takes the display
settings too.

**The viewer** is split into a `ViewState` -- table, spatial index, filter,
settings, and a field of view in, an image out -- and a matplotlib front end, so
another toolkit replaces only the second half.  Three things make it usable on
millions of localizations:

* it renders at *display* resolution, so the pixel size follows from the window
  and the cost of a render is bounded by the canvas rather than by the zoom;
* a gesture only moves the axes limits and matplotlib rescales the image it
  already has; a timer re-renders once the gesture has been still for 150 ms.
  Interaction stays smooth even where a render takes 130 ms;
* the `SpatialIndex` answers what is on screen, so a zoomed view never walks the
  whole table.  Queries are padded by `roi_sigma * max_sigma` so a blob whose
  centre is off screen still contributes its tail -- the index may return extra
  candidates but never loses one, and that is what the tests check.

The filter controls cover only the fields that get used (localization precision,
z or PSF size, `logl_rel`, frame) rather than every column; a slider at either
end means "no bound", so opening a file never silently drops localizations.

Threading is by contiguous row bands: a localization contributes to whichever
bands its kernel reaches and each band writes only its own rows, so there are no
shared pixels, no locks and no per-thread image copies.

## Grouping

A single emitter is localized in many consecutive frames, so a raw table counts
one blink many times.  `group.py` links them (SMAP's `connectsingle2c.c`, greedy
frame-to-frame, first match inside a **box** of half-width `dx` wins, up to `dt`
dark frames) and combines each run into one row using the per-column rules from
`Grouper.m`.

* **The linking bug is fixed.**  The original tested `list[thisentry] > 0`
  before the bounds check, so with the last localization already assigned it
  read -- and could then write -- one past the end of the array.  `Grouper.m`
  carries two workarounds for the unassigned first/last entries ("FIX
  connectsingle doesnt assign last loc. Fix later!"); with the tests in the
  right order every localization gets a group and the workarounds are gone.
* **Blocks are linked separately** rather than by SMAP's trick of zeroing the
  frame at each boundary, which leaves the array no longer sorted by frame and
  lets a group start in one block and search into the next.
* **The error of a summed quantity adds in quadrature.**  This is the rule that
  keeps shot noise consistent: with `e_i = sqrt(N_i)`, `sqrt(sum(e_i^2))` is
  `sqrt(sum(N_i))` exactly -- the shot noise of the summed photons.  SMAP
  applies its `*err -> 1/sqrt(sum(1/e^2))` rule to `photons_err` and
  `background_err` instead, which understates a four-member group by six times.
  Background is summed alongside photons, so the pair stays consistent: a group
  is one measurement over several exposures, with total signal and total
  background.  (Note this makes a grouped `background` n times the per-frame
  level; a filter on it means something different than on ungrouped data.)
* **`logl` is dropped, `logl_rel` is kept.**  A group's raw log-likelihood is a
  sum over its members' fits, so no reduction of it is meaningful across groups
  of different size; the per-pixel `logl_rel`, which is what the filter uses,
  takes the maximum as in SMAP.
* **Each coordinate is weighted by its own error** -- `x_nm` by `x_err_nm`,
  `y_nm` by `y_err_nm`, `z_nm` by `z_err_nm` (`WEIGHT_FOR`).  SMAP weights all
  three with the pooled `locprecnm` because that is the one weight it computes,
  and flags it as a shortcut.  The pooled weight is the correct inverse-variance
  estimate only when the errors are equal, and under astigmatism they are not:
  on the clathrin data `x_err/y_err` runs from 0.42 to 3.15 across
  localizations, since that divergence is exactly what encodes z.  The
  difference in the grouped position is a median 0.11 nm laterally but reaches
  14-21 nm, and in z a median 0.41 nm reaching 487 nm -- far above the
  localization precision, so it is not cosmetic.

Grouped and ungrouped are **two tables kept side by side**, each with its own
filter and spatial index (`LocSet`), because a precision cut that makes sense
for single localizations is the wrong one for groups, and because rebuilding the
filter caches on every switch would throw away the work.  The viewer's switch is
then a single assignment; only the first switch pays for the linking.

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
* the render kernel's second moment is `sigma^2 + 1/12`, the variance the pixel
  itself contributes -- a point-sampled kernel would give `sigma^2`
* the C++ renderer matches a numpy reference to float32 precision in all three
  modes, coloured and not
* the whole renderer runs on real data: 844k localizations from the clathrin
  dataset, filtered to 575k, render to structure that looks like clathrin pits,
  and zooming to a 150 nm field resolves the individual localizations as
  ellipses of their own precision
* the spatial index is checked against brute force on random rectangles: it may
  return extra candidates, never lose one
* grouping on the clathrin dataset: 844k -> 377k rows, photon sum conserved
  exactly, `frame` is the group start, and the precision of every multi-member
  group is better than its best member (median 7.0 -> 4.5 nm)
* the summed-error rule is checked against the shot-noise identity: for
  `e_i = sqrt(N_i)` the grouped error is `sqrt(sum(N_i))` to float precision

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

### Rendering

M1 Pro, 5M localizations into 2000x2000:

| | full field (100 nm/px) | zoomed 10x (10 nm/px) |
|---|---|---|
| histogram | 49 ms | 23 ms |
| fixed sigma | 156 ms | 40 ms |
| sigma from precision | 133 ms | 42 ms |
| the same, coloured | 314 ms | 48 ms |
| the same, one thread | 898 ms | 44 ms |

Filtering 5M localizations on four fields: 18 ms to build, 5 ms for the AND,
3 ms to compact.  Moving one slider costs 14 ms (one field recomputed, AND,
compact) plus 31 ms to render the 0.7M that survive -- live filtering at ~20 fps.

Grouping 844k localizations (dx = 50 nm, dt = 1): **0.31 s**, dominated by the
sequential linking.  The population median precision *rises* (17.5 -> 21.0 nm)
even though every group improves, because bright long-lived emitters collapse
into few rows and dim singles stay -- the per-group comparison above is the one
that means something.

End to end through the viewer, 5M localizations filtered to 3.7M, into a
1000x1000 canvas: **126 ms** for the full field, **3.9 ms** zoomed 10x, **0.4 ms**
zoomed 100x (select + render).  Without the index the zoomed cases cost ~42 ms,
because every thread walked all five million localizations to find the ones on
screen.

The spatial index over 5M localizations builds in 100 ms and answers in
0.02-0.12 ms.  It is a `lexsort` on 16-bit row and column keys: numpy radix-sorts
narrow integer keys and falls back to a comparison sort for wide ones, so the
obvious `argsort` on a combined cell id costs 1.5 s instead of 60 ms.

Rendering at *display* resolution rather than a fixed nm/pixel is what bounds
this: zoomed out the kernels collapse onto the `mingausspix` floor, zoomed in
there are fewer localizations in view.  The zoomed numbers are floored by every
thread scanning all 5M localizations to find the ~50k in view; a spatial index
would remove that.

## Two windows

The image has a window to itself and the controls a second, small one.  The
image window is then free to be as large as the screen allows, and since the
rendered pixel size follows the canvas rather than the zoom, enlarging it
renders a *finer* image rather than scaling one up -- on the clathrin file, the
same view went from 60.3 to 35.5 nm per rendered pixel simply by not sharing
the window with a control strip.  The controls keep their size instead of
stretching across a wide window, and their rows are laid out in inches with the
window sized to fit, so the spacing does not change with the number of filter
fields or colour choices.  Closing the image closes both; closing the controls
leaves the image alone.

A side effect worth having: the text boxes are no longer on the image canvas,
so a click on the image cannot reach them at all (see below).

The image window may be any shape.  The axes is never shrunk to the proportions
of the view -- `set_aspect("equal", adjustable="datalim")` -- so it is the view
that widens to the proportions of the window: a wide window shows more x rather
than growing bands beside the image.  One pixel size still serves both axes,
which is not negotiable; an SMLM image with non-square pixels would be a lie.
`FieldOfView.fit` already worked this way (it sizes a render from the canvas
and picks the pixel size from the larger span), so the two agree by
construction: on a 1300x600 window the render is 1248x558 pixels of 46.5 nm,
covering 58.0 x 26.0 um.  Empty area beside square data stays black, as part of
the image, rather than grey figure background.

Filter bounds now have defaults -- `loc_precision_nm` below 25 nm, `logl_rel`
above -1.5, `z_nm` within +-500 -- because those are right almost every time
and an unfiltered first view is not what anyone wants to look at.  Every one of
them throws localizations away, so each is written into its box: the earlier
rule that opening a file must never *silently* drop localizations still holds,
it is the silence that mattered, not the absence of a default.  They apply once
per table (the grouped table gets its own), a bound already set is never
overridden, and a bound the user clears stays cleared.  Only fields with a
fixed unit get one, since 25 nm means nothing in pixels.

## Interaction latency

Two things made a drag feel like it started a second late, neither of them the
render (35 ms for 211k localizations at this canvas size):

**Every `TextBox` redraws the whole figure when a click lands anywhere else**,
to take its cursor away -- with `canvas.draw()`, which blocks.  A window with a
dozen boxes therefore paid a dozen full draws before the first motion event of
a drag was delivered: 212 ms measured here, close to a second on a retina
canvas.  `viewer._patched_text_box` swaps `draw` for `draw_idle` around
`stop_typing`, collapsing all of them into the one draw that was going to
happen anyway; the press went from 212 ms to 0.9 ms.  It is swapped in around
the call rather than reimplementing `stop_typing`, whose body reaches into
private state that changes between matplotlib versions.

**A queued render could land in the middle of the gesture.**  The press now
stops the timer.

With those gone, deferring the render to the end of a gesture is no longer
worth it when the render is cheap: `_refresh` renders inside the gesture when
the last one took under `LIVE_RENDER_SECONDS` (0.12 s), throttled by the
render's own cost so a gesture never queues up renders it cannot keep up with,
and falls back to moving the image and deferring when it is slower -- a slow
render inside a gesture would make it lurch, which is worse than a blank edge.
Panning and zooming therefore fill in what they expose as they go.  `fit`
keeps the centre and is idempotent, so re-rendering repeatedly inside a drag
cannot make the view creep: measured drift over a twelve-event drag is 7e-12
nm.

## Online analysis

Fitting an acquisition while the microscope writes it and watching the image
build up (`smapfit.live`, `smapfit.io.watch`).

**The bridge between the fitter and the viewer is a queue in one process, not
the HDF5 file.**  The file is still written, and it is still the result -- but
using it as the transport would mean SWMR mode (every column declared before
the first block knows what they are, `refresh()` on every read) and re-reading
and decompressing a tail that is already sitting in memory.  The engine already
hands finished blocks to a `sink`; the live path adds a second sink, which is a
queue.  A separate process only earns that cost when the fitter and the viewer
are on different machines, and then the change is additive: same
`ViewState.append`, a different feeder.

**The fitter runs in a thread, and the GIL is not in the way.**  Detection,
fitting, rendering and grouping all release it in the C++ (`csrc/`), and frame
reading already overlaps through `prefetch`, so the window stays responsive
while a block is being fitted.  Only the viewer's own thread touches the view,
on a timer, so nothing in the state needs a lock: a block of localizations is
handed over and never looked at again by the thread that made it.

**An update must not change a setting.**  This is structural, not a discipline:
appending *extends* the table, each filter's cached mask and the spatial index,
and re-derives none of them.  What that took:

* `Localizations.extend` grows each column in a buffer with 50% headroom and
  exposes an exact-length view, so 500 updates over an hour do not each copy
  the whole table.
* `LocFilter.append` evaluates the bounds already set on the new rows only.  A
  bound therefore keeps meaning exactly what was typed; nothing re-quantiles
  the larger table.  A hand-set mask has no rule to extend it, so its new rows
  are excluded until it is set again -- visible, rather than silently included.
* `GrowingIndex` gives each appended block its own `SpatialIndex` and answers a
  query as the union.  Segments are merged occasionally (the tail when there
  are too many of them, everything when the tail reaches the size of the head),
  so the merge work is O(N log N) over an acquisition rather than a 0.3 s
  argsort per update.  Segments deliberately do not share a grid: a query is
  answered against each segment's own origin, so a block landing outside the
  area seen so far needs no rebuild.
* The index takes an `extent` -- the camera's field of view -- so `bounds`, and
  with it the viewer's full view, is the same before the first localization and
  after the last.  The image does not rescale under the user.
* The median localization precision, which caps the rendering sigma, is
  recomputed once the table has grown by a fifth rather than on every block.

**Colouring by a field is a re-render, and its range is explicit.**  Each
localization contributes its own colour, so the LUT is applied per localization
and baked into the accumulation -- unlike contrast and gamma, which re-display
what is already there.  The range is never taken from the data at render time,
which `render_locs` would otherwise do: the visible set changes with every zoom,
every filter and, online, every block, so an implicit range would give the same
z a different colour each time.  Choosing a field therefore fills the range in
from that field's 1--99% quantiles and leaves it there; a range typed for one
field is not carried to the next, since a z slab means nothing for a frame
number.  The LUT follows the choice (`hot` for intensity, `turbo` for a field),
because a coded field cannot be read off a brightness ramp.

**Grouping is the one thing that cannot be extended.**  It links localizations
across frames, so an append marks the grouped table stale and it is rebuilt the
next time it is asked for, keeping the bounds set on it.  The title says
"(stale)" rather than the view quietly lagging.

**Knowing when the acquisition ended.**  Micro-Manager does not say, and the
declared frame count is what was planned.  So: read what is there, wait, look
again, and stop when nothing new has appeared for `timeout` seconds.  The
newest page of the file being written is held back until another follows it
(its pixels may still be arriving) and read on the way out, once nothing can
follow it.  Every read is retried rather than trusted -- a file caught
mid-write raises instead of returning half an image -- and the file list is
re-globbed on each poll, so the `_1.ome.tif` that appears when the current file
fills up is picked up without reopening anything.

**The engine is flushed on a timer as well as when its ROI buffer fills.**  The
buffer holds 15000 ROIs; a sparse sample would take minutes to fill it, and
nothing would appear in the meantime.

Checked against the offline path on 300 frames of the astigmatic dataset,
replayed frame by frame into a growing two-file series: same 28724
localizations, every column identical.

## Drift correction

Estimation is not ours: COMET (vendored under `externaltools/Comet`) is called
as a library. `drift.py` is the seam between it and the rest of the package,
and it is deliberately thin -- an array in, a per-frame drift table out.

* **Estimate on a subset, correct everything.** The estimator wants clean
  localizations; a badly fitted one contributes noise to the overlap cost and
  nothing else. So `estimate_drift` takes the same `select` a render does. But
  the correction is a coordinate change and reaches every row, filtered out or
  not -- otherwise changing a filter afterwards would mean re-running it.
* **COMET's file I/O is not wrapped.** It can load ThunderSTORM CSV and write
  molecule sets; we hand it a `(N, 4)` array and read the drift back, so the
  corrected file stays an ordinary smapfit HDF5 that the viewer opens unchanged.
  The drift curve rides along in a `/drift` group, which readers ignore.
* **The drift table is indexed by frame.** COMET interpolates onto
  `arange(0, max_frame + 1)`, so applying it is one fancy-index per coordinate.
  `min_max_frames` is passed explicitly because the subset the estimate is made
  from need not reach the last frame of the table it is applied to.
* **`peak_*` is left alone.** Those are detection positions -- the record of
  what was found where -- not measurements of the structure.
* **Units survive.** COMET works in nm; a table in pixels is converted for the
  estimate and the drift divided back by the pixel size when applied.
* **A new file, not an in-place edit.** `_driftc.hdf5` beside the input: drift
  correction is a hypothesis (segmentation width, maximum drift) that is worth
  redoing with other parameters, and the raw fit must stay untouched to redo it.

### Where the time goes

Measured on the clathrin dataset (410 k localizations after filtering, 92 time
windows, 314 M neighbour pairs):

| stage | time |
|---|---|
| segmentation | 0.03 s |
| pair-count estimate (the safety check) | 0.6 s |
| KD-tree neighbour search | 5 s |
| **optimization** | **97 s** (was 309 s) |
| applying the drift and writing the file | ~5 s |

The optimization is everything: L-BFGS-B evaluates cost and gradient a few
hundred times, and each evaluation is one `exp()` per pair over all 314 M of
them -- ~0.6 s, compute-bound at ~2 ns/pair across the cores.  So there are only
two levers: the work per evaluation (the pair count) and the number of
evaluations.  Two changes to the vendored COMET, both opt-in from `drift.py`
(`CPU_WRAPPER` and `LBFGSB_OPTIONS` in `comet/core/drift_optimizer.py` keep
upstream's values for anyone importing COMET directly):

* **`cpu_wrapper_chunked_fast`** -- the original converts its arguments on every
  call, and the pair indices arrive as int32 and are copied to int64: 5 GB of
  allocation and copying per cost evaluation, for nothing, since the compiled
  kernel specialises on whatever dtypes it gets.  0.73 -> 0.61 s per call, same
  answer to 5e-10.
* **`ftol = 1e-7`** instead of COMET's `1e3 * eps` (~2e-13).  423 -> 154 cost
  evaluations.  The drift changes by 0.12 nm rms (2.2 nm on one window out of
  92) and the rendered image by 0.03% -- against a localization precision of
  10-20 nm, that tolerance was buying nothing.

Together: 309 s -> 97 s (and 75 s with the approximate kernel below).  Things that were tried and did **not** help, so that
they are not tried again:

* float32 coordinates, int32 indices, and sorting the localizations spatially so
  that paired ones sit near each other in memory: no change at all.  The kernel
  is compute-bound on `exp`, not waiting on memory, which is also why the
  gradient scatter-add is not worth restructuring.
* a chunked NumPy version, to get NumPy's SIMD `exp`: **44x slower** -- it is
  single-threaded and materialises the per-pair temporaries.
* dropping same-segment pairs (they add a constant and exactly zero gradient):
  correct, but only 6.8% of the pairs.
* a smaller `max_drift_nm`.  Pairs grow with r^3 only for uniform data; this is
  clustered, and going 300 -> 150 nm removes just 21% of the pairs.
* the PyTorch/MPS backend on Apple silicon: the index tensors alone are 5 GB on
  the device and a single cost evaluation did not finish in ten minutes.

### Two more that were tried

* **A looser tolerance for the coarse sigma steps** -- the idea being that they
  only have to land the estimate in the right basin.  **Refuted.**  Coarse
  `ftol=1e-4` is 2.3x faster but moves the drift by 4.8 nm rms and 52 nm on
  individual windows, and a *tighter* final step does not repair it: a local
  optimizer at small sigma cannot undo a wrong basin.  Loosening `gtol` instead
  traces the same curve (48 evaluations, 2.4 nm rms), so it is not an artefact
  of `ftol` being a relative criterion.  Anything below ~140 evaluations costs
  tens of nm.  The knob survives as `optimizer_ftol_coarse`, defaulting to off.
* **An approximate kernel** (`cpu_wrapper_chunked_approx`).  Kept: **1.35x**,
  100.8 -> 74.8 s end to end, drift unchanged (0.12 nm rms against the reference,
  the same as the exact kernel) and image sharpness identical to four digits.
  `exp` is 54% of the kernel on ARM, so it is tabulated on the integer part and
  Taylor-expanded on the fraction; pairs beyond 6 sigma are skipped, and the
  four divisions per pair become one hoisted reciprocal.  A 4 sigma cutoff was
  measured too: no faster, and two orders of magnitude less accurate.

### The noise floor, and how to measure it

Deviation from a previous run is a bad yardstick without knowing how much a
drift estimate varies by chance.  Split the filtered localizations into two
disjoint halves, estimate from each, and their difference is pure noise by
construction; since noise goes as 1/sqrt(N), the full estimate's noise is half
the difference between the halves:

| estimate | halves differ by | noise of the full estimate |
|---|---|---|
| ungrouped | 3.25 nm rms | **1.6 nm** |
| grouped | 6.19 nm rms | **3.1 nm** |

So a change that moves the drift by less than ~2 nm rms is not measurable, which
is what makes `ftol=1e-7` (0.12 nm) and the approximate kernel (0.12 nm) safe,
and what makes the coarse-tolerance experiment (4.8 nm) a real degradation.

### Grouping before estimating

`DriftSettings.group` estimates from grouped localizations -- one entry per
blink instead of one per frame.  On the clathrin dataset 410 k -> 146 k
localizations, and **314 M -> 21 M pairs**: the whole correction takes 5 s
instead of 1:17.  It beats the 2.8^2 that the localization count alone suggests,
because it removes exactly the dense very-short-range pairs, which is where the
pair count grows fastest.

Grouping links in x and y only and averages z, which is right: one localization
per frame is one molecule (multi-emitter fits are removed by the log-likelihood
cut), and a molecule does not move in z between consecutive frames, so averaging
z only makes it more precise.

**Filter before estimating, and include z.**  Without a z cut the two estimates
disagreed by 8.2 nm rms and the ungrouped z drift range came out 113 nm against
the grouped 75 nm -- which looked like grouping damaging the axial estimate, and
was not.  `z_err_nm` has a median of 41 nm and a 95th percentile of 108 nm: the
out-of-focus tail carries the axial estimate away.  With
`logl_rel > -2`, `loc_precision_nm < 15`, `|z_nm| < 300` (323 k of 844 k
localizations) the ranges agree (69 vs 75 nm) and the two curves lie on top of
each other: **median |difference| 0.9 nm in x, 0.6 in y, 1.4 in z**.

What is left is that grouping doubles the noise of the estimate (3.1 nm against
1.6 nm, split-half) because there are fewer, though better, localizations per
window.  On the image that costs a few percent: of the improvement the full
estimate makes over no correction, grouping recovers 93% in the x-y projection
and 82% in x-z.  So the full estimate stays the default and grouping is the
switch for a first look, for large datasets, and wherever 35x matters more than
a few percent of resolution.

### A failed window, and quality control on the CPU

The ungrouped estimate on this dataset contains one **bad time window** (49,
centred on frame 24 819): the drift swings 46 nm in y and 34 nm in z and comes
straight back, over 2500 frames, where the grouped estimate -- same data, more
precise positions -- moves 6 nm.  A stage does not do that.  It is what inflates
the grouped-vs-ungrouped rms (4.3 nm in y) over its median (0.6 nm).

COMET's answer to this is quality control: per segment, compare the overlap the
fitted drift achieves against the overlap with no drift correction at all, and
NaN the segments that fail, so the spline bridges them.  It was wired into the
`cuda_qc` and `torch_qc` backends only; the CPU path raised "not implemented".
It now works there (`cpu_wrapper_chunked_qc`, `DriftSettings.quality_control`,
`--qc`), one extra pass over the pairs, about 4 s of the 150 s run.

**The criterion, however, does not work as written, and the replacement only
half works.**  Both findings are worth keeping:

* Upstream flags a segment when ``q_obs < q_null + std(q_null)``, where the
  standard deviation is taken *across* segments.  That spread is not an
  uncertainty on any one segment: on a bleaching sample the pair count per
  window falls from 14 M to 2 M, so the null varies five times more between
  segments than the correction improves any of them.  It flagged **90 of 92**
  windows here; with 18 left the spline bridged gaps of thousands of frames and
  the estimated drift range blew up from 53/75/69 nm to 491/267/591 nm.
* Normalising instead by the null itself -- the **lift**,
  ``q_obs / q_null - 1`` (`flag_flawed_segments_by_lift`, the default here) --
  removes that density dependence, and at sigma = 20-30 nm it does single out
  window 49 as the only segment whose fitted drift makes its own overlap
  *worse* than no correction.  But the lift is a function of sigma, and at the
  sigma this run actually finishes at (13 nm) window 49 comes out at +0.06 and
  is not flagged.  So the flag is real but marginal, and it catches this failure
  only at some sigmas.

The reason neither works well is that overlap-versus-null is not the quantity
one wants.  Window 49 has few pairs (2.7 M against 14 M early in the movie), so
its cost landscape is flat: the estimate is badly *constrained*, while its
overlap is not obviously bad.  What that asks for is the **curvature of the cost
in each segment's three parameters**, which is an uncertainty per window in nm,
accumulated in the same single pass and needing a threshold in nm rather than a
tuned dimensionless one.  That is the thing to implement next; it would also
give the drift curve error bars, which nothing currently does.

### Redundant cross-correlation

`smapfit.rcc` is the other classical estimator -- bin the acquisition into time
windows, render each, and measure how far one has moved against another by
cross-correlation -- ported from SMAP's `finddriftfeature.m` /
`finddisplacementZ2.m`.  It shares no code and almost no assumptions with the
overlap estimator, which is the point: it is the second opinion on a drift
curve, and it takes seconds.

*Redundant* means every pair of windows is correlated, not just consecutive
ones, giving T(T-1)/2 measurements of `d_k - d_l` for T unknowns; the drift is
the robust least-squares solution of that overdetermined system, so a
correlation peak found in the wrong place is outvoted instead of accumulated.

Compared on the same footing -- same filter, all estimating from **grouped**
localizations, each judged on the half of the molecules it did not see:

| estimator | time | noise x/y/z (nm) | out-of-sample x-y | x-z |
|---|---|---|---|---|
| **COMET spline, 2000 fr knots (default)** | 3.6 s | 0.69 / 0.61 / 0.73 | **4.4249e-5** | **3.7097e-5** |
| COMET spline, 2300 fr knots (= RCC window) | 2.0 s | 0.63 / 0.58 / 0.61 | 4.4178e-5 | 3.7068e-5 |
| COMET per-window | 1.7 s | 1.25 / 1.70 / 1.45 | 4.3983e-5 | 3.6563e-5 |
| RCC, 20 windows | 7.4 s | 0.61 / 0.90 / 0.86 | 4.3843e-5 | 3.6580e-5 |
| no correction | | | 3.5268e-5 | 3.1192e-5 |

At matched smoothing -- spline knots at the RCC window spacing -- the two agree
to **1.2 / 1.7 / 1.6 nm rms**, at the level of their own noise, from methods
sharing no code and no assumptions.  The spline is the better of the two on both
counts here, but not by much, and RCC's axial noise is close to it.

Grouping matters as much to RCC as it does to COMET, so it is on by default
there too.  Comparing a grouped estimator against an ungrouped one, which is
what the first version of this table did, mostly measures the grouping.

Two things were needed to get there, both of which cost a factor of ten or more
in accuracy when wrong:

* **The axial correlation must be tiled in both x and y.**  SMAP's
  `finddisplacementZ` slices only in x, each slice spanning the whole field;
  each z profile then pools hundreds of structures into a slab whose correlation
  peak has nothing to lock onto.  Measured: pairwise shifts *uncorrelated* with
  COMET's (-0.09) with x-slices, +0.96 with 200 nm square tiles.  The tiled
  version, `finddisplacementZ2`, is the one to follow.
* **The lateral drift must be taken out before the axial pass**, or a tile holds
  different structure in different time windows and the z profiles being
  correlated are not of the same thing.

`exclude_zero_lag` leaves the zero-lag sample out of the axial peak search and
fit, on the grounds that anything sitting at the same z in both windows --
a molecule on in both, localizations pinned to a z plane by the fit -- pulls the
estimate towards zero drift.  Measured here it makes no difference (0.86 nm
axial noise excluded against 0.70 nm kept, identical agreement with COMET),
because grouping already collapses each molecule to one entry and this data has
no pinning: no 5 nm z bin holds more than three times the median count.  It
stays on as insurance for data that does; a z histogram is how to check.

Other departures: the correlation peak is located by a quadratic fit rather than
a 2D Gaussian (the peak is not Gaussian, and this needs no starting values); the
overdetermined system is solved by iteratively reweighted least squares rather
than MATLAB's robust `nlinfit`, since `d_k - d_l` is linear in the unknowns and
needs no nonlinear fit; per-frame interpolation is PCHIP, which does not
overshoot between windows.

### Drift as a spline in time

`DriftSettings.spline` (`--spline`) replaces the free vector per time window
with a cubic B-spline in time, fitted directly: the optimizer's variables are
the spline's coefficients, the cost and gradient are COMET's over the same
pairs, and the chain rule to the coefficients is one matrix product.  There is
**no segmentation** -- each localization sits at its own frame -- and no
interpolation afterwards, because the fitted curve *is* the per-frame drift.
B-splines are variation-diminishing, so it cannot overshoot the way an
interpolating spline through noisy window estimates does.

Drift is slow creep with the occasional small jump, so a coefficient every 2000
frames (26 of them here, against 276 free window parameters) describes it, and
each coefficient is constrained by every localization near it in time rather
than by one window's worth.

Judged on halves that are **disjoint in molecules** (see below), and by
correcting each half with the *other* half's drift so that no estimate is graded
on localizations it was fitted to:

| estimator | time | noise x/y/z (nm) | out-of-sample x-y | x-z |
|---|---|---|---|---|
| windows, ungrouped | 34 s | 1.87 / 3.73 / 5.19 | 94.0% | 100.0% |
| windows, grouped | 2 s | 1.25 / 1.70 / 1.45 | 100.0% | 90.0% |
| spline, ungrouped | 31 s | 1.62 / 3.39 / 5.03 | 94.1% | 100.7% |
| **spline, grouped** | **2 s** | **0.93 / 0.81 / 0.99** | **102.9%** | 95.9% |

Knot spacing on the full data: 200 frames gives 1.7-2.0 nm noise, 2000 gives
0.6-0.7 nm, and the out-of-sample sharpness is flat from 1000 to 4000 -- no sign
of over-smoothing, which is what a slowly drifting stage should look like.
`spline_penalty` adds a second-difference penalty on top if a dataset needs it.

### Two ways to be wrong about which estimate is better

Both of these produced a wrong conclusion earlier in this work, so they are
written down:

* **Split by molecule, not by localization.**  A blink is several localizations
  of one molecule in consecutive frames.  Splitting localizations at random puts
  one molecule in both halves, the halves share their errors, and the noise
  comes out too low -- it said the ungrouped estimate had a 1.6 nm noise floor
  against the grouped one's 3.1 nm.  Splitting whole molecules reverses it:
  1.9/3.7/5.2 nm ungrouped against 1.3/1.7/1.5 nm grouped.  Grouping does not
  throw information away; treating a molecule's repeats as independent
  measurements adds correlated noise, and over-weights bright molecules.
* **Grade out of sample.**  Rendering the same localizations the drift was
  fitted to rewards the estimator that fitted their noise.  In sample the
  ungrouped estimate looked sharper (3.9761e-5 against 3.9056e-5); correcting
  each half with the other half's drift reverses that laterally.

The remaining honest gap is axial: the ungrouped estimate is still slightly
better in x-z (100% against 95.9%), even though its z drift is by far the
noisiest.  Worth understanding rather than explaining away.

### Two passes: grouped, then a small radius

`DriftSettings.two_stage` (`--two-stage`) runs the grouped estimate first, takes
that drift out, and then runs an ungrouped estimate over a **30 nm** search
radius.  The pair count is what costs time, and once the coarse pass has removed
the drift the fine pass only has to look a few nm: 314 M pairs become 42 M.

| | time | vs the single pass (x/y/z) | sharpness x-y | x-z |
|---|---|---|---|---|
| single pass, ungrouped | 153 s | reference | 3.9761e-5 | 3.3623e-5 |
| grouped only | 6 s | 1.8 / 4.3 / 4.7 nm | 3.9056e-5 | 3.2388e-5 |
| grouped + fine, r=100 nm | 91 s | 1.5 / 6.6 / 6.3 nm | 3.9465e-5 | 3.3769e-5 |
| grouped + fine, r=50 nm | 26 s | 1.2 / 4.0 / 3.5 nm | 3.9583e-5 | 3.3626e-5 |
| **grouped + fine, r=30 nm** | **18 s** | 0.8 / 3.6 / 3.0 nm | **3.9743e-5** | 3.3513e-5 |

**8.7x faster for 99.8% of the lateral and 98.4% of the axial improvement**, and
theremaining difference is at the noise floor (1.6 nm).  Note that the earlier estimate
of what this would gain -- 2.5-3.5x -- was too pessimistic: it assumed the fine
pass would need as many cost evaluations as the single pass, and it does not,
because it starts a few nm from the answer instead of a hundred.

It is also **more robust**, which was not the point of it.  The fine pass is
bounded by its own radius, so it cannot make a large excursion: over the window
where the single pass fails (frames 24 000-26 500) the single-pass estimate
swings 46 nm in y, and the two-stage one 7.9 nm.  A small search radius is a
speed knob and a regularizer at once.

The radii are tuned on one dataset.  `two_stage_radius_nm` sets the fine pass's
radius and, from it, its sigma schedule (radius/3 down to radius/5); it must
stay well above the residual the grouped pass leaves, which was a few nm here.

What is left is a real trade, not a free win, so it stays a parameter:
`max_locs_per_segment` caps the localizations per time window and the pair count
falls roughly quadratically -- 2000 per window is 10x faster (33 s) and moves
the drift by 2.5 nm rms, 1000 per window is 15x faster and moves it by 6 nm.
Worth it for a first look at a large dataset, not for the final result.

## Open questions

* Fitted x sits ~0.24 px from the peak-finder position, and the sign flips with
  the mirror (-0.22 unmirrored, +0.26 mirrored) while y stays at -0.09.  A
  round trip (render the model, fit it back) is exact, so this is a property of
  the data/calibration, not of the code.  Worth comparing against SMAP's own
  output on the same dataset.
* Camera offset: taken from the image metadata, but an iXon does not report one,
  so the settings file is used with a warning.  Should that be a hard error?
* RCC's slice width, z bin width and number of time windows are SMAP's values,
  which were set by hand rather than optimized; there is likely room there.
* `group()` numbers groups from 1, a leftover from the MATLAB original, so the
  row of the grouped table is `group_index - 1`.  Worth rebasing to 0 at some
  point, together with the C++ `connect`.
* Why the ungrouped estimate still wins slightly in the axial projection when
  its z drift is three times noisier than the grouped one's.
* A configuration layer that merges file metadata with user settings once, up
  front, so downstream code never deals with partial metadata.
* A filter bound cannot be set on a live table before the first block has
  arrived, because the column does not exist yet: `LocFilter.set` raises and
  the viewer restores the box.  Holding such a bound as pending and applying it
  on the first append would be friendlier.
* Online drift correction, and an online estimate of the labelling density:
  both want the whole table, which is what an appended view already has.

## Environment

`python3` on this machine is x86_64 conda under Rosetta.  Use the native arm64
venv:

    .venv/bin/python -m pytest tests/
    .venv/bin/python setup.py build_ext --build-lib src

The extension builds universal2.  `setup.py` contains a workaround for a
Command Line Tools install with incomplete libc++ headers; it is gated on a
test-compile and stays silent on a healthy toolchain.
