# pyCOMET

**pyCOMET** (Cost-function Optimized Maximal overlap drift EsTimation) is a method and Python package for high-precision drift correction in Single Molecule Localization Microscopy (SMLM) data.

---

## Background

Fluorescence imaging techniques such as **STORM, PALM, DNA-PAINT, and MINFLUX** push spatial resolution in biological specimens below 10 nm. These approaches rely on recording the spatial coordinates of fluorophores over time. Inevitably, **mechanical drift** caused by thermal expansion or mechanical instabilities shifts fluorophore positions relative to the optical system, which distorts localization data and reduces image resolution.

Existing approaches to drift correction include:

- **Hardware solutions** (drift-stable microscope designs)  
- **Fiducial marker tracking** (accurate, but requires suitable markers in the focal plane)  
- **Correlation-based methods** (marker-free, operate directly on localization coordinates, but limited in time resolution and speed)

Each of these methods carries limitations in terms of practicality, accuracy, or computational efficiency.

---

## Concept

COMET introduces an **optimization-based approach** that requires no fiducial markers and is capable of detecting drift with high accuracy and temporal resolution.

The core assumptions are:

1. A unique time-varying displacement vector **D(t)** exists that maximally overlaps all localizations when subtracted.  
2. Given sufficient sampling at each time point, this displacement corresponds to the **true sample drift**.  

COMET formulates drift correction as the **minimization of a cost function** that quantifies the overlap of localizations. This cost function is derived from a Gaussian mixture model:

- Each localization is represented as a Gaussian distribution.  
- The cost reflects the summed overlap across all localization pairs.  
- Minimization aligns localizations over time by solving for D(t).

---

## Implementation

Direct evaluation of the cost function over all localization pairs is computationally prohibitive. COMET addresses this by:

- Restricting calculations to pairs within a maximum drift distance (**Dmax**)  
- Running cost function and gradient evaluation on the **GPU (CUDA)** or CPU fallback  
- Using an iterative **L-BFGS-B optimizer**, combined with a scale-reduction schedule of the Gaussian parameter (σ) to avoid local minima and converge robustly  

This iterative strategy solves drift at progressively finer scales, yielding both robustness and precision.

---

## Features

- **3D drift correction** without fiducials  
- GPU-accelerated for high-throughput datasets  
- Configurable temporal resolution via **segmentation modes** (frames, localizations, or windows)  
- Flexible interpolation of segment-wise drift to per-frame drift  
- Integration via both **Python API** and **CLI tool**

---

## Summary

COMET provides a **fast, marker-free, optimization-based solution** to drift correction in SMLM datasets. By working directly on localization coordinates and leveraging GPU acceleration, it achieves higher accuracy and temporal resolution compared to correlation-based methods, while maintaining computational efficiency.

