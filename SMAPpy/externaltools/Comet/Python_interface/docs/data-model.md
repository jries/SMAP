# Data model (in‑memory)

## Localization array
- Shape: `(N, 4)` → columns `[x_nm, y_nm, z_nm, frame]`
- Units: nm; `frame` is integer
- Source: `drift_optimizer.py`

## Segments
- Frames are mapped to segment IDs `0..S−1` before optimization.
- Invalid/sampled‑out localizations are masked out.
- Source: `segmenter.py`, `drift_optimizer.py`.

## Per‑segment drift
- Vector `μ ∈ ℝ^{S×3}` (x,y,z per segment).
- Interpolated to **per‑frame drift** at the end.
- Source: `interpolation.py`, `drift_optimizer.py`.

> **Frames vs segments**: The corrected localization array returned by
`comet_run_kd(..., return_corrected_locs=True)` uses **segment IDs** in the last column (not original frames).
If you need per‑frame drift, request the drift curve and map frames → drift yourself.