# Changelog

All notable changes to the COMET Python package are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

COMET is released in two lanes:

| Lane      | `requires-python` | Purpose                                                       |
| --------- | ----------------- | ------------------------------------------------------------- |
| **1.1.x** | `>=3.9`           | Active development.                                            |
| **1.0.x** | `>=3.6`           | Frozen lane so pip does not block Python 3.6/3.7 installs.     |

`pip install py-comet` resolves to the right lane automatically: pip skips any
release whose `requires-python` excludes the running interpreter. Python 3.6/3.7
are not a supported configuration — the lane exists so those users are not
refused at install time.

## [Unreleased]

## [1.1.0] - 2026-08-17

First release published to PyPI.

### Added

- Automatic backend detection (`comet.best_backend()`, `comet.available_backends()`,
  `comet.describe_backends()`). `comet_self_test` and the `comet` CLI now pick the
  fastest available backend instead of assuming CUDA.
- A public API surface on the `comet` package itself: `from comet import comet_run_kd`,
  the IO helpers, and `comet.__version__`.
- `comet_self_test` gained `--mode cpu`, `--size {quick,full}` and `--verbose`, reports
  the detected environment, and returns a meaningful exit code.
- `interactive=` parameter on `comet_run_kd` for the large-pair-count confirmation
  prompt, off by default.
- A test suite of 167 tests that runs without a GPU (`pytest --pyargs comet.tests`).
  GPU-dependent tests are marked and skip automatically.
- `torch`, `test`, `docs` and `dev` extras.

### Changed

- **The `cpu` backend is now compiled with numba, making it 380–500x faster.**
  It was a pure Python loop over neighbour pairs. numba was already a hard
  dependency — it is what the CUDA kernels are built on — so this costs nothing
  extra to install and needs no GPU.

  End-to-end pipeline, unchanged accuracy:

  | Localizations | Before  | After  | Speed-up |
  | ------------- | ------- | ------ | -------- |
  | 1 500         | 3.4 s   | 0.23 s | 15x      |
  | 3 200         | 19.7 s  | 0.04 s | 460x     |
  | 6 000         | 88.1 s  | 0.17 s | 505x     |

  The full-size `comet_self_test` on CPU alone now finishes in 13 s rather than
  several minutes, and the project's own test suite dropped from 121 s to 22 s.

  Above 200 000 pairs a parallel kernel is used, which is a further ~3x on an
  8-core machine; the threshold is where the two were measured to cross over.
  The plain-Python implementation is retained as `_cost_and_gradient_reference`
  and is the specification the compiled kernels and GPU backends are tested
  against.

- The CPU path no longer allocates an `(n_pairs,)` float64 scratch buffer it
  never needed — 800 MB at 100 million pairs.

- **`best_backend()` no longer selects torch on an Apple MPS machine.** With the
  compiled CPU kernel in place, MPS is ~6x *slower* for this workload across the
  realistic range — 3.4 vs 23.8 ms at 496 k pairs, 19.7 vs 108.6 ms at 4 M,
  57.0 vs 349.5 ms at 12.5 M. The MPS GPU itself is fine (2.8 TFLOP/s on
  matmul); eager-mode torch just materialises every intermediate across all
  pairs, so at 4 M pairs each `(P, 3)` temporary is ~48 MB and the evaluation is
  bound by memory traffic. The compiled kernel keeps each pair in registers.
  `torch_accelerated()` is therefore CUDA-only now.

  This is a speed decision, not an accuracy one: MPS is float32-only and its
  per-evaluation gradients differ by up to ~1e-5 at 3 M pairs, but the converged
  drift agrees with the CPU backend to 0.0013 nm.

  `mode="torch"` still selects MPS explicitly; it is only the automatic choice
  that changed. `describe_backends()` now says why.

- **`comet_run_kd(mode=...)` now defaults to `None`, meaning auto-detect**, instead of
  `"cuda"`. On a CUDA machine the selection is unchanged; on every other machine the
  previous default raised `CudaSupportError`, so no working code changes behaviour.
  This is what made the documented Python API example fail for most users.
- **`requires-python` is now `>=3.9`.** Python 3.6/3.7 users are served by the 1.0.x
  lane; see the table above. The previous `">=3.6,<=3.14"` was wrong at both ends —
  it excluded Python 3.14.1, and 3.6 was never actually installable.
- Dependency caps `numpy<2.0` and `matplotlib<3.8` removed. Verified against
  numpy 2.5, scipy 1.18, pandas 3.0 and matplotlib 3.11.
- Per-iteration optimizer progress is no longer printed unconditionally; it is gated
  behind the existing `display` flag.
- `plot_q_with_baseline()` no longer computes the flawed-segment indices; that moved to
  the new `flag_flawed_segments()` so quality control can run headless.

### Removed

- `scikit-learn` and `lmfit` dependencies, which were declared but never imported.

### Fixed

- `import comet` failed on any Python without tkinter (common on Linux and in
  containers), because four modules imported `tkinter` at module level. The file
  dialogs now load on demand.
- `comet_self_test` failed on every machine without an NVIDIA GPU.
- `comet --format h5 -o out.h5` ignored `--output` and opened a save dialog.
- `comet.utilities.post_analysis_utilities` imported the undeclared dependency `sympy`,
  which also broke `comet.batch.mfx_batch`.
- Restored `analyse_folder_of_h5_drift_summary_files()` and
  `create_nice_plot_from_x_and_y_drift()`, removed accidentally in 378b5f0.
- `cuda_qc` mode raised `NameError` when torch was not installed, because `qc_utils`
  was imported inside the torch guard.
- The `disp` L-BFGS-B option, removed in newer SciPy, raised `OptimizeWarning` on every
  optimizer call.
- A segmentation parameter too coarse for the dataset (more frames per window than the
  movie has frames) failed with `ValueError: 'x' must contain at least 2 elements` from
  inside SciPy. COMET now reports what went wrong and how to fix it.
- The README documented `--initial_sigma_nm` default 600, `--target_sigma_nm` default 1
  and `--max_drift_nm` default None; the actual CLI defaults are `max_drift_nm/3`, 30
  and 100.
- The citation link in the README had the DOI URL duplicated inside it, so it did not
  resolve.
- **Molecule-set HDF5 output was wrong in three ways.** `save_dataset_as_ms_h5()` wrote
  the first coordinate column to `Y_POS_PIXELS` and the second to `X_POS_PIXELS`, so a
  save/load round trip returned x and y swapped. It also divided z by the lateral pixel
  size instead of the axial one, scaling z by `pixelsize_nm / pixelsize_z_nm`. And
  `(N, 2)` input was written with `Z_POS_PIXELS = 1`, which loaded back as a spurious
  z offset of one pixel instead of zero.
- `load_normal_molecule_set()` checked the legacy general pixel-size key before the
  newer axis-specific pair, so a file carrying both scaled z by the general value.
  The axis-specific keys now take precedence, a lone `xy_pixel_size_um` falls back for
  z, and pixel sizes stored as length-1 arrays are handled.
- **`save_corrected_locs=True` saved uncorrected localizations.** It wrote
  `sorted_dataset`, the segmented and possibly downsampled copy taken *before* the
  drift was subtracted, whose last column holds segment ids rather than frame numbers.
  The saved molecule set therefore contained raw positions with segment ids in
  `FRAME_NUMBER`. It now writes the corrected array with real frame numbers.
- The molecule-set pixel size was hardcoded to 160 nm on that path. It is now the new
  `pixelsize_nm` / `pixelsize_z_nm` arguments to `comet_run_kd()`, exposed on the CLI as
  `--pixelsize_nm` / `--pixelsize_z_nm`. Molecule sets store positions in pixel units,
  so this must match the acquisition for the coordinates to mean anything outside COMET.
- `comet_run_kd`'s docstring described the returned `corrected_locs` last column as
  `segment_id`; it is the frame number.
- **`interpolation_method="catmull-rom"` returned integer drift.** The result array was
  allocated with `np.zeros_like(frame_range)`, and `comet_run_kd` passes an integer
  frame range, so every drift value was truncated to whole nanometres — a 3.7 nm drift
  became 1 nm. Frames outside the interior spans were also left at exactly zero,
  reading as "no drift" at the start and end of the acquisition; they are now clamped
  to the nearest in-range estimate.
- **The CUDA backend reported an inflated cost on datasets above ~100 M pairs.**
  `d_val_sum` accumulates on the device across chunks and is never reset, but the host
  loop added each intermediate copy as well, weighting chunk *k* by `n_chunks - k + 1`.
  L-BFGS-B was therefore given a cost inconsistent with its gradient, which can show up
  as line-search failures on large datasets. Single-chunk runs were unaffected.
- `comet --format csv` with an `.h5` input passed the HDF5 path to `pd.read_csv`. It now
  writes the corrected array directly.
- `estimate_pairs()` and `pair_indices_lex_floor_asymmetric()` shifted the caller's
  coordinate array to the origin in place.
- `temporal_refined_drift()` in `comet.core.experimental` unpacked two values from
  `pair_indices_kdtree()`, which returns three.
- `save_correction_details=True` always stored `gt_drift=None`, discarding the ground
  truth the caller supplied.
- `correct_and_save_thunderstorm_csv()` raised `IndexError` when the drift table was
  shorter than the CSV's frame range; it now explains the mismatch.
- The CLI's input-extension check was case-sensitive, rejecting `DATA.CSV`.
- `pair_indices_kdtree()` returned bare lists on `MemoryError`, so callers doing
  `.size` or `.astype` on the result failed with a confusing `AttributeError`.

### Changed (from the dev branch)

- CLI `--max_drift_nm` now defaults to 300 nm, matching the library default (was 100).

## [1.0.2]

Compatibility lane for Python 3.6/3.7. Same source as 1.1.0, with dependency floors
lowered to the last releases supporting those interpreters and the conditional
`dataclasses` backport added.

[Unreleased]: https://github.com/gpufit/Comet/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/gpufit/Comet/releases/tag/v1.1.0
[1.0.2]: https://github.com/gpufit/Comet/releases/tag/v1.0.2
