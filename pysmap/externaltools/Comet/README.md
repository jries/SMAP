![image >](Python_interface/resources/comet_logo_small.png)

**Cost-function Optimized Maximal Overlap Drift EsTimation**
[Preprint available](https://doi.org/10.64898/2026.03.27.714864)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/gpufit/Comet/blob/master/Colab_notebooks/COMET.ipynb)

## Overview

**COMET** is a fast, GPU-accelerated software package for drift correction in single-molecule localization microscopy (SMLM) datasets. It achieves high spatial and temporal resolution by maximizing spatiotemporal overlap across frames using a cost-function optimization approach.

---

## How to Use

### 1. Try Online via COMET Web Tool

Visit our web platform at [smlm.tools](https://www.smlm.tools), upload your dataset, and get results directly without any installation required.

#### Prepare your file:

* Format: CSV (ThunderSTORM-compatible)
* Required headers: `"frame"`, `"x [nm]"`, `"y [nm]"`, and optionally `"z [nm]"`
* Column headers must match exactly (quotes included)
* Extra columns are allowed

[ThunderSTORM reference](https://zitmen.github.io/thunderstorm/)

#### Upload and Run:

1. Upload your file on [smlm.tools](https://www.smlm.tools)
2. Choose:

   * Segmentation method (e.g. segment by number of localizations per time window)
   * Segmentation parameter (e.g. 500 locs per time window)
   * The maximum drift expected in nm
3. Click **Run**

#### Tips:

* ✅ Check *"Keep file for later"* to re-analyze with different settings
* ✅ Check *"Spline Interpolation"* to get a smooth result per frame
* ✅ Check *"Dynamic downsampling"* if you experience memory errors on very large datasets
* ⏳ Busy? If queue times are high, use the Python/Colab version locally.

---

### 2. Run Locally (Python Package & CLI)

**Requirements**

* Python 3.9 or newer (3.6/3.7 are still supported — see [Legacy Python](#legacy-python-36--37))
* Optionally a CUDA-capable GPU compatible with Numba CUDA, for full acceleration
* NumPy, SciPy, Matplotlib, Pandas, h5py and Numba are installed automatically

**A GPU is optional.** COMET runs on three backends and picks the fastest one
available automatically:

| Backend | Needs                              | Notes                                        |
| ------- | ---------------------------------- | -------------------------------------------- |
| `cuda`  | NVIDIA GPU + CUDA driver           | Fastest; the numba-cuda kernels              |
| `torch` | `pip install "py-comet[torch]"`    | Uses CUDA or Apple MPS, or falls back to CPU |
| `cpu`   | nothing                            | numba-compiled; no GPU needed                |

#### Installation

```bash
pip install py-comet
```

For the PyTorch backend (useful on machines without an NVIDIA GPU, including
Apple Silicon):

```bash
pip install "py-comet[torch]"
```

To test the installation:

```bash
comet_self_test
```

This simulates a dataset with a known drift, recovers it, and reports which
backends it found. It needs no GPU and takes a few seconds. Add `--plot` to see
the recovered drift, or `--mode cpu` to force a specific backend.

<details>
<summary>Installing from source instead</summary>

```bash
git clone https://github.com/gpufit/Comet
cd Comet/Python_interface
pip install -e ".[dev]"
```

This creates an "editable" install so changes to the code take effect
immediately. See [CONTRIBUTING.md](CONTRIBUTING.md) for the development
workflow.

</details>

<details id="legacy-python-36--37">
<summary>Legacy Python (3.6 / 3.7)</summary>

Python 3.6 and 3.7 are **not supported, but deliberately not blocked**. The aim
is that pip does not refuse the install, so you can try COMET on an older
interpreter if that is what you have available.

`pip install py-comet` resolves to the frozen `1.0.x` release line on those
interpreters, whose dependency floors are low enough for them. You do not need a
git URL or a different package name.

Caveats:

* Upgrade pip first (`python -m pip install --upgrade pip`). Versions older than
  pip 9 ignore the `requires-python` metadata that makes this work, and would
  try to install a release that needs 3.9+.
* Getting `numba` and `llvmlite` installed on these interpreters is fiddly and
  platform-dependent — expect to match wheels by hand.
* CI only covers this path on a best-effort basis, so treat it as "try it and
  see". Use Python 3.9+ if you have the choice.

The `1.0.x` lane receives backported fixes only; new features go to `1.1+`.

</details>

#### CLI Usage

Once installed, you have a `comet` command available in your environment:

```bash
comet \
  --input       your_data.csv \
  --output      corrected.csv \
  --format      csv \
  --segmentation_mode 2 \
  --segmentation_var 60
```

To see all options:

```bash
comet --help
```

#### Key CLI Parameters

* `--input` (string, **required**): Path to your input file (`.csv` or `.h5`).
* `--output` (string, **required**): Where to save the corrected output.
* `--format` (csv|h5, default=csv): Output file format.
* `--pixelsize_nm` (float, default=160): Camera pixel size in nm, recorded in `--format h5` output. Molecule sets store positions in pixel units, so set this to match your acquisition.
* `--pixelsize_z_nm` (float, default: same as `--pixelsize_nm`): Axial pixel size for `--format h5` output.
* `--segmentation_mode` (0|1|2, default=2):

  * `0`: Fixed number of time windows (`--segmentation_var` = number of segments)
  * `1`: Fixed number of localizations per window (`--segmentation_var` = locs per segment)
  * `2`: Fixed number of frames per window (`--segmentation_var` = frames per segment)
* `--segmentation_var` (int, **required**): Value associated with your chosen mode.
* `--initial_sigma_nm` (float, default: `max_drift_nm / 3`): Initial Gaussian sigma (nm) for overlap optimization.
* `--target_sigma_nm` (float, default=30): Target sigma (nm) at which the algorithm stops refining.
* `--max_drift_nm` (float, default=300): Maximum expected drift (nm). Also the neighbour-search radius.
* `--boxcar_width` (int, default=1): Width of the moving-average filter applied between iterations.
* `--interpolation` (cubic|catmull-rom, default=cubic): Interpolation method for per-frame drift curves.
* `--max_locs_per_segment` (int, default=None): Cap on localizations per time window, to bound memory and runtime.
* `--mode` (cuda|torch|cpu|cuda_qc|torch_qc, default: fastest available): Compute backend.
* `--display`: Print progress and show diagnostic plots.

You can omit any optional parameters to use their default values.

---

### 3. Python API

If you prefer to call COMET directly in Python:

```python
from comet import (
    comet_run_kd,
    load_thunderstorm_csv,
    save_dataset_as_thunderstorm_csv,
)

# Load your CSV -> (N, 4) array of [x_nm, y_nm, z_nm, frame]
dataset = load_thunderstorm_csv("your_data.csv")

# Run drift correction
drift, corrected = comet_run_kd(
    dataset,
    segmentation_mode=2,
    segmentation_var=60,
    initial_sigma_nm=100,
    target_sigma_nm=10,
    max_drift_nm=300,
    boxcar_width=1,
    interpolation_method="cubic",
    return_corrected_locs=True,
)

# Save as CSV
save_dataset_as_thunderstorm_csv(corrected, "corrected.csv")
```

`comet_run_kd` returns `drift` as an `(F, 4)` array of
`[dx_nm, dy_nm, dz_nm, frame]`, one row per frame.

**Choosing a backend.** The default is `mode="cuda"`. To let COMET choose, or to
check what is available:

```python
import comet

print(comet.describe_backends())          # what this machine supports
drift = comet_run_kd(dataset, mode=comet.best_backend(), ...)
```

> **Note:** `comet_run_kd` modifies the array you pass in — it subtracts the
> estimated drift from `dataset` in place. Pass `dataset.copy()` if you need to
> keep the original.

---

### 4. Google Colab Notebook

Click the badge above to launch the interactive version in your browser with no setup required.

---


## Documentation

This repository also ships developer documentation built with [MkDocs](https://www.mkdocs.org/) and the [Material theme](https://squidfunk.github.io/mkdocs-material/).  
It includes usage guides, background, and an auto-generated API reference.

### Build locally

Install the documentation extras:

```bash
pip install "py-comet[docs]"
```

Then from within the `Python_interface` folder build and serve the docs:

```bash
mkdocs serve
```

Open your browser at http://127.0.0.1:8000

## Segmentation Modes

COMET segments data before estimating drift. Choose from:

| Mode | Description                             | Parameter                 |
| ---- | --------------------------------------- | ------------------------- |
| 0    | Fixed number of time windows            | Number of segments        |
| 1    | Fixed number of localizations/window    | Localizations per segment |
| 2    | Fixed number of frames/window (default) | Frames per segment        |

---

## Experimental Features
Some features are experimental and may change in future releases:
### GPU acceleration for Mac (MPS) and AMD GPUs (ROCm)
GPU acceleration for Mac (MPS) and AMD GPUs (ROCm) is under development.
Using pytorch a first running version of COMET on these platforms is possible, but performance may vary. 
Using the same (cuda capable) GPU initial tests showed that the numba cuda implementation is currently 
at least 2x faster than the pytorch implementation. Anyhow, to enable usage of COMET on MPS and ROCm platforms
we included a pytorch based implementation. First successful tests were done on Apple Silicon (M2). 
Feedback from users with AMD GPUs is welcome! 

#### Installing PyTorch for MPS
To install PyTorch with MPS support on macOS, use the following command:

```bash
pip install torch torchvision torchaudio
```
To test if the installation was successful and comet can run using torch try the following command:

```bash
comet_self_test --plot --mode torch
```
then simply call COMET as usual with the mode input parameter specified to 'torch'.  

## Citation

> If you use COMET in your research, please cite this [preprint](https://doi.org/10.64898/2026.03.27.714864).
> Machine-readable citation metadata is in [CITATION.cff](CITATION.cff).

---

## Contact

For questions or contributions, feel free to open an issue or reach out on GitHub.
