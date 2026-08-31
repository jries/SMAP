# I/O helpers

## CSV (ThunderSTORM)
`load_thunderstorm_csv(...)`
- Returns `(N,3)` if no `"z [nm]"` column → `[x, y, frame]`
- Returns `(N,4)` if z is present

**2D tip**: For 3D pipeline, add zero‑z: `np.c_[x, y, np.zeros_like(x), frame]`.

## Molecule set (HDF5)
- Load: `load_normal_molecule_set(...)`
- Save: `save_dataset_as_ms_h5(...)`

## Save correction details
`save_drift_correction_details(...)`
- HDF5 groups with drift curves + parameters
- Some functions open Tk dialogs when `filename=None` (non‑headless friendly)