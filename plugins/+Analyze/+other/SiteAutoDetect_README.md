# SiteAutoDetect — SMAP Plugin for Automatic CCP Site Detection

## What it does

`SiteAutoDetect.m` is a SMAP plugin that automatically detects clathrin-coated pit (CCP) sites from SMLM localizations and adds them directly to the SMAP Site Explorer (SE). No manual clicking required.

## Algorithm

```
All localizations (xnm, ynm)
        ↓
Render 2D density map  (10 nm/pixel, accumarray)
        ↓
Gaussian blur  (blur_sigma_nm)
        ↓
Normalise to 99th percentile  (better contrast than max-normalisation)
        ↓
imregionalmax()  →  keep only local maxima above threshold
        ↓
Greedy min-distance suppression  (min_dist_nm)
        ↓
Centroid refinement  (recenter each peak on the centroid of locs within se_siteroi)
        ↓
Count locs per site  (within se_siteroi radius)
        ↓
Apply min-locs filter  (exclude OR mark with Use=false)
        ↓
Assign sites to a grid of cells  (grid size = se_cellfov from SE Settings)
        ↓
Add sites + cells to SE
```

Works for both **2D and 3D** data (only xnm/ynm are used for detection; z is set to 0).

## Parameters

| Parameter | Default | Description |
|---|---|---|
| Threshold | 0.02 | Fraction of the 99th-percentile density map value. Lower = more (fainter) sites detected. |
| Min site distance (nm) | 150 | Minimum distance between two detected peaks. Prevents double-detection of the same CCP. |
| Blur sigma (nm) | 40 | Gaussian smoothing of the density map before peak finding. Should roughly match the apparent size of one CCP spot. |
| Min locs per site | 10 | Sites with fewer localizations within `se_siteroi` are filtered. |
| Locs filter action | Exclude | What to do with sites below min locs: **Exclude** = don't add to SE; **Mark (Use=false)** = add with `-` in the site list. |
| Histogram bins | 30 | Number of bins for the locs-per-site histogram shown in Preview. |
| Clear existing sites | unchecked | If checked, removes all existing sites (and cells) from the SE before adding new ones. |

`se_cellfov` and `se_siteroi` are read automatically from the SE Settings panel (not shown in the plugin GUI).

## Buttons

- **Preview** — runs the full pipeline and shows a figure with:
  - Left panel: density map with detected sites overlaid (cyan = pass, yellow × = below locs filter)
  - Right panel: histogram of locs per site with a red dashed line at `min_locs`
  - Does **not** modify the SE.
- **Detect Sites** — same pipeline, writes results to the SE.
- **Run** (SMAP's built-in button) — same as Detect Sites.

## Tuning tips

Start with Preview. Typical good starting values for CCP data:
- `blur_sigma = 40 nm` (keep low — CCPs are small and dense)
- `threshold = 0.02–0.05`
- `min_dist = 150 nm` (CCPs are ~150–200 nm apart minimum)
- `min_locs = 10–30` depending on labelling efficiency

Use the locs histogram (Preview) to decide where to set `min_locs` — look for the natural separation between real sites and noise peaks.

## Files

| File | Location | Purpose |
|---|---|---|
| `SiteAutoDetect.m` | `plugins/+Analyze/+other/` on the `ccp-autodetect` branch | Plugin class |
| `SiteAutoDetect_README.md` | same folder | This document (install + usage instructions) |

`SiteAutoDetect.m` is available at:
**https://github.com/AndreuBoixPages/SMAP/tree/ccp-autodetect/plugins/%2BAnalyze/%2Bother**

## Installation

This plugin has been submitted to the official SMAP repository as **PR #34**:
**https://github.com/jries/SMAP/pull/34**

Once merged, just pull the latest SMAP and the plugin appears automatically — no manual steps needed:

```bash
cd <your SMAP folder>
git pull
```

Then in MATLAB: `clear classes` (or restart MATLAB), and the plugin will appear under **Analyze → other → Auto-detect CCP Sites**.

> **Before the PR is merged:** download `SiteAutoDetect.m` from
> `https://github.com/AndreuBoixPages/SMAP/tree/ccp-autodetect/plugins/%2BAnalyze/%2Bother`
> and copy it manually to `<SMAP_root>/plugins/+Analyze/+other/`. No changes to `plugin.m` needed —
> SMAP auto-discovers plugins at startup.

## After detection — verifying it worked

```matlab
se = g.locData.SE;
fprintf('Sites detected: %d\n', length(se.sites));
fprintf('First site at (%.0f, %.0f) nm\n', se.sites(1).pos(1), se.sites(1).pos(2));
```

## Known limitations / future work

- Cell grid uses a fixed tile size (`se_cellfov`). For very non-uniform datasets, cells at the edge of the FOV may be only partially filled.
- Centroid refinement uses a circular ROI of radius `se_siteroi`. If two CCPs are closer than this, the centroid may be pulled between them. Consider reducing `se_siteroi` in the SE Settings if CCPs are very dense.
- Currently uses only 2D positions for detection. A future version could incorporate z-position to detect 3D clusters in thick samples.
- The min-locs filter counts all localizations within the ROI regardless of z. For 3D data, adding a z-range filter would improve specificity.
