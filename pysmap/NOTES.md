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
