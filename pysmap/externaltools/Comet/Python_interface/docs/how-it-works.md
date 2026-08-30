# How it works

COMET formulates drift correction as an **optimization problem**:  
find the time-dependent displacement vector D(t) that maximizes the overlap of all molecular localizations in a dataset.

---

## Cost function

Each localization Ai is shifted by its segment-wise drift vector:

Ai'(D(ti)) = Ai - D(ti)

The pairwise distance is:

d_ij(D(t)) = || Ai'(D(ti)) - Aj'(D(tj)) ||

The COMET cost function is defined as:

fC({Ai}, D(t)) = - Σ_ij exp( - d_ij^2(D(t)) / (2 * sigma^2) )

- Localizations are modeled as Gaussian kernels with width sigma.  
- Maximizing overlap = minimizing distances between corrected localizations.  
- The gradient of fC can be calculated analytically, enabling efficient optimization.

---

## Optimization strategy

- **Segmentation:** The dataset is split into temporal windows (by frames, by localizations, or fixed number of windows).  
- **Pairs:** Spatially close pairs are identified (KD-tree).  
- **Backend:** Cost function and gradients are computed on GPU (CUDA) or CPU fallback.  
- **Minimizer:** L-BFGS-B is used to optimize D(t).  
- **Schedule:** sigma is decreased iteratively (coarse → fine), ensuring robust convergence to the global minimum.  
- **Smoothing:** Optional temporal boxcar averaging regularizes the drift trajectory.

---

## Interpolation

The optimized per-segment drift vectors mu (S × 3) are interpolated to per-frame drift curves:

- **Cubic spline** (`scipy.interpolate.CubicSpline`)  
- **Catmull–Rom spline** (requires ≥4 points)  

This yields a continuous 3D drift trajectory over all frames.

---

## Output

- **Per-frame drift curve**: array of shape (F × 4), with columns `[dx, dy, dz, frame]`  
- **Corrected localizations**: optional, with drift subtracted (last column stores segment IDs, not original frame numbers)  

---

## Key advantages

- Fiducial-free (requires only localization coordinates)  
- Sensitive to both slow drift and fast transients  
- Scales efficiently with GPU acceleration  
- Flexible: API functions and CLI command available