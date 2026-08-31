# How COMET works (step by step)

This page explains the full processing flow in **plain text**, using the actual function names in the codebase.
It is intended as a conceptual "flowchart in words" that you can follow alongside the source.

---

## Big picture

Goal: estimate a time-dependent drift D(t) that, when subtracted from all localizations, makes them maximally overlap.

High-level stages:

1) **Prepare data** -> ensure array shape (N, 4) == [x_nm, y_nm, z_nm, frame]  
2) **Segment in time** -> map frames to segment IDs 0..S-1  
3) **Find close pairs** -> KD-tree radius graph within drift_max_nm  
4) **Optimize per-segment drift** -> minimize cost function using L-BFGS-B  
5) **Smooth and schedule** -> optional boxcar smoothing, reduce sigma over iterations  
6) **Interpolate** -> per-segment drift mu(S x 3) -> per-frame drift (F x 3)  
7) **Optionally correct locs** -> subtract drift, write outputs (CSV/HDF5)

All of this is orchestrated by `comet.core.drift_optimizer.comet_run_kd(...)`.

---

## Data model and entry

- Input array: `(N, 4)` with columns `[x_nm, y_nm, z_nm, frame]` (nm units, frame is int).
- Function: `comet.core.drift_optimizer.comet_run_kd(dataset, segmentation_mode, segmentation_var, ...)`
- Typical preprocessing:
  - If dataset has no z column, add zeros: `np.c_[x, y, np.zeros_like(x), frame]`.
  - Keep frames as recorded; segments are derived later.

**Why this first?** The optimizer and KD-tree expect consistent units and a fixed column order.

---

## 1) Temporal segmentation

Function: `comet.core.segmenter.segmentation_wrapper(loc_frames, segmentation_var, segmentation_mode, ...)`

- **mode 0**: `segment_by_num_windows` (choose S directly)  
- **mode 1**: `segment_by_num_locs_per_window` (accumulate frames until >= X locs)  
- **mode 2** (default): `segment_by_frame_windows` (fixed window size in frames)

**Outputs (SegmentationResult):**
- `loc_segments (N,)`: segment ID per localization (−1 when unused)
- `loc_valid (N,)`: mask of kept localizations (sampling may drop excess locs)
- `center_frames (S,)`: representative frame per segment
- `n_segments`: S
- `out_dict`: counts/metadata

**Why segmentation?**  In principle the intrinsic temporal resolution of any data-based dirft correction in SMLM is fundamentally limited to a single frame, so ideally a drift estimate per frame is obtained. However, the typical SMLM dataset frame does not contain enough information to extract a reasonable drift estimate on a per frame basis, making it necessary to bin frames together to obtain a sane result --> temporal segmentation. 
The optimizer then estimates one drift vector per segment. Good temporal bins balance time resolution and per-bin information. A good starting point is usually a few hundred localizations per segment.

---

## 2) Pair search (KD-tree radius graph)

Function: `comet.core.pair_indices.pair_indices_kdtree(coordinates, distance)`

**Inputs:**
- `coordinates = locs_nm[:, :3]` (x,y,z in nm) for the **kept** localizations
- `distance = drift_max_nm` (usually 3 * initial_sigma_nm, unless overridden)

**Outputs:**
- `idx_i, idx_j` (int32 arrays) listing pairs (row-wise)

**How it works:**
- Build a KD-tree over 3D points, query pairs within the given radius.
- If memory becomes an issue, the radius can be reduced internally or the dataset can be downsampled manually (optinal functionality included in the segmentation step) to prevent crash.

**Why now?** The cost function quantifies the overlap of **nearby** pairs and finding these pairs is computationally expensive. Computing pairs once and reusing them across optimization steps makes the entire algorithm fast. Even if in this initial step not all localizations belonging together are paired up (e.g. due to huge drift) and a lot of pairs are found that should not be overlapped the information of the correct pairs is enough to find the true drift. 

---

## 3) Backend selection and memory layout

Files: `cuda_wrapper.py`, (optinally but not recommended: `cpu_wrapper.py`)  both wrapped via `drift_optimizer.py`

**Setup:**
- Allocate arrays for drift parameters `mu` (shape `(S, 3)`) and intermediate buffers on device (GPU) or host (CPU).
- Transfer `idx_i, idx_j`, coordinates, segments to device if using GPU.

**Why this step?** Proper memory layout and a single transfer of static data (pairs, coordinates) avoids repeated copies and keeps the inner optimization loop fast.

---

## 4) Cost function and gradient

Ai'(D(ti)) = Ai - D(ti)
d_ij = || Ai'(D(ti)) - Aj'(D(tj)) ||
fC = - sum_{i,j} exp( - d_ij^2 / (2 * sigma^2) )
grad = analytic derivative of fC w.r.t. per-segment drift mu


Implementation:
- GPU kernels: `cost_function_full_3d_chunked` in `cuda_wrapper.py`
- CPU kernels: `cost_function_full_3d_chunked_cpu` in `cpu_wrapper.py`
- Both accumulate the **negative sum** (objective) and its **analytic gradient** using atomic adds per segment.

---

## 5) Optimizer and schedule

Function: `comet.core.drift_optimizer.optimize_3d_chunked_better_moving_avg_kd(...)`

**Loop outline:**
1. Start with larger `sigma` (coarse scale) to find global structure.
2. Run several L-BFGS-B steps on flattened `mu` (shape `S*3`), subject to box constraints:
   - bounds: `+- drift_max_nm * drift_max_bound_factor` (to keep updates physical and 'help' the optimizer)
3. Between steps, apply optional **boxcar smoothing** (temporal regularization on `mu`).
4. If improvement slows, reduce `sigma` by a factor (typ., `/ 1.5`), and continue from the current `mu`.
5. Stop near `target_sigma_nm` or when objective no longer improves.

**Why this order?**
- Coarse-to-fine `sigma` avoids local minima in the cost function landscape early on.
- Bounds prevent extreme, non-physical solutions and helps the optimizer to converge fast.
- Smoothing encourages temporal consistency without overfitting noise.

---

## 6) Interpolation to per-frame drift

Function: `comet.core.interpolation.interpolate_drift(center_frames, drift_est, frame_range, method='cubic'|'catmull-rom')`

- Inputs: `center_frames (S,)` and `drift_est mu (S, 3)`
- Outputs: per-frame drift `(F, 3)` for frames `0..max_frame`
- Methods:
  - `cubic` -> `scipy.interpolate.CubicSpline`
  - `catmull-rom` -> requires at least 4 segment points

`comet_run_kd` then builds `frame_interp = np.arange(0, max_frame+1)` and returns:
- `drift_interp_with_frames` as `(F, 4)` with columns `[dx, dy, dz, frame]`

**Why interpolate?** The optimizer estimates per-segment drift; users often want **per-frame** drift curves. 
This of course assumes a smooth drift between the time windows and might not be the ideal thing to do in all cases!

---

## 7) Optional: subtract drift and save

- If `return_corrected_locs=True`, `comet_run_kd` also returns corrected localizations:
  - Note: the last column holds **segment IDs** (not original frames).
- Save helpers in `comet.core.io_utils`:
  - `save_drift_correction_details(...)` -> HDF5 groups (drift curves + parameters)
  - `save_dataset_as_ms_h5(...)` -> molecule set HDF5
  - CSV output for curves or corrected locs is also supported in CLI

**Tip:** For headless use, always pass explicit filenames to avoid GUI dialogs.

---

## Pseudocode (end-to-end)

```
def comet_run_kd(dataset, segmentation_mode, segmentation_var, **kw):
    # 0) validate input shape (N,4)
    locs = ensure_shape_n4(dataset)  # [x_nm, y_nm, z_nm, frame]

    # 1) temporal segmentation
    seg = segmentation_wrapper(
        loc_frames=locs[:, 3].astype(int),
        segmentation_var=segmentation_var,
        segmentation_mode=segmentation_mode,
        max_locs_per_segment=kw.get("max_locs_per_segment"),
    )
    keep = seg.loc_valid.astype(bool)
    locs_kept = locs[keep]
    seg_ids = seg.loc_segments[keep]

    # 2) KD-tree pairs within radius = drift_max_nm
    coords = locs_kept[:, :3]
    drift_max_nm = kw.get("max_drift", 3 * kw.get("initial_sigma_nm", 100))
    idx_i, idx_j = pair_indices_kdtree(coords, drift_max_nm)

    # 3) choose backend, allocate, transfer once if GPU
    backend = "cuda" or "cpu" (or "torch", currently experimental)
    state = backend_prepare(backend, coords, seg_ids, idx_i, idx_j)

    # 4/5) optimize mu with schedule on sigma
    mu = optimize_3d_chunked_better_moving_avg_kd(
        n_segments=seg.n_segments,
        locs_nm=coords,
        seg_ids=seg_ids,
        pairs=(idx_i, idx_j),
        initial_sigma_nm=kw.get("initial_sigma_nm", 100),
        target_sigma_nm=kw.get("target_sigma_nm", 1),
        boxcar_width=kw.get("boxcar_width", 1),
        backend_state=state,
        bounds_scale=kw.get("drift_max_bound_factor", 1.0),
    )

    # 6) interpolate to per-frame
    max_frame = int(locs[:, 3].max())
    frames = np.arange(0, max_frame + 1, dtype=int)
    drift_interp = interpolate_drift(
        center_frames=seg.center_frames,
        drift_est=mu,          # (S,3)
        frame_range=frames,
        method=kw.get("interpolation_method", "cubic"),
    )
    drift_with_frames = np.c_[drift_interp, frames]

    # 7) optional corrected locs (store segment IDs in last column)
    if kw.get("return_corrected_locs", False):
        corrected = locs.copy()
        corrected[keep, :3] = corrected[keep, :3] - mu[seg_ids]
        corrected[:, 3] = seg.loc_segments  # segment ids
        return drift_with_frames, corrected

    return drift_with_frames

````

Why this order works

Segmentation before pairs: defines temporal bins so pairs can be reused while mu is updated.

KD-tree once: pair graph is static, so the expensive neighbor search is not inside the optimizer inner loop.

Coarse-to-fine sigma: encourages global alignment first, then fine-scale refinement, improving convergence and stability.

Interpolation last: consumers typically need per-frame drift, not per-segment parameters.

### Parameter checklist

- segmentation_mode and segmentation_var: pick so each segment has enough localizations. Start near 200 locs/segment.
- initial_sigma_nm -> how strong a pair of localizations feels each other given a distance they are apart, large values make the cost function landscape smooth, will be refined stepwise.
- target_sigma_nm -> refine to 1..5 nm for final precision (data dependent, not equivalent to spatial resolution of drift estimate, e.g. in simulations a spatial resolution of the drift estimate < 1nm could be reached with a target sigma value of 10 nm).
- max_drift -> default: 3 * initial_sigma_nm; increase if expected drift is larger.
- boxcar_width -> temporal smoothing of mu; 0 or 1 for minimal smoothing, larger for noisy data.
- interpolation_method -> "cubic" is robust; "catmull-rom" needs >= 4 segments.
- mode -> default is "cuda", set cpu on systems without CUDA (rather for testing purposes than day-to-day use). In the future setting it to torch will allow usage of non-nvidia GPUs but that's currently an experimental feature.

### Common pitfalls
- Wrong input shape: use (N,4) exactly; for 2D CSVs, insert zero z.
- Too few locs per segment: optimizer may stall or become unstable. Increase segment size.
- Overly small max_drift: true pairs might be missed or optimizer hits bounds; increase radius.
- Forgetting that corrected_locs last col == segment IDs: do not treat it as original frame index.
- Dialogs in headless runs: pass file paths to save functions to avoid GUI prompts.

### Where to look in code
- Orchestrator: comet/core/drift_optimizer.py (comet_run_kd, optimizer entry)
- Segmentation: comet/core/segmenter.py
- KD-tree pairs: comet/core/pair_indices.py
- CUDA kernels: comet/core/cuda_wrapper.py
- CPU kernels: comet/core/cpu_wrapper.py
- Interpolation: comet/core/interpolation.py
- I/O helpers: comet/core/io_utils.py
- CLI wrapper: cli_interface.py