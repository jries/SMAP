# fit3Dcspline
fit3Dcspline is a GPU based 3D single molecule fitter for arbitrary, experimental 
point spread functions (PSF). The fitting speeds achieves more than 10^5 fits/s on a 
consumer graphic card GTX 1070. The implmentation of the fitting algorithm is based
on maximum likelihood estimation and employs Levenberg-Marquardt optimization routine, 
which reaches theoretical minimum uncertainty. Both the EMCCD and sCMOS noise model 
are included. The softare package also includes tools to robustly model beads based
experimental PSFs of different modality and correct for depth induce aberrations. 

# Requirements
  - Matlab R2016a or newer  
    - Curve Fitting Toolbox
    - Optimization Toolbox

The GPU fitter requires:
  
  - Microsoft Windows 7 or newer, 64-bit
  - CUDA capable graphics card, minimum Compute Capability 3.0
  - CUDA 8 compatible graphics driver (for GeForce products 378.66 or later)

The CPU version runs on macOS and Microsoft Windows 7 or newer, 64-bit
  
# How to run
Example code for 3D PSF calibration based on beads on coverglass is avalible in file **example_3D_fit.m**. Example code for depth dependent PSF calibration based on beads in gel is avalible in file **example_depth_aberration.m**. The required 3D image stacks for
the demo code can be found in the folder by following [this link](https://oc.embl.de/index.php/s/L9nZGqFzfIIzps8). A full instruction guide can be found 
in **User_guide_Ries.pdf**.

# Contact
For any questions / comments about this software, please contact [Ries Lab](https://www.embl.de/research/units/cbb/ries/index.html).

# Copyright and Software License
Copyright (c) 2017 Ries Lab, European Molecular Biology Laboratory, Heidelberg. 

fit3Dcspline also includes OME Bio-Formats package for reading and converting biological
file formats in folder **bfmatlab** which comes with a separate copyright. 

The fit3Dcspline is licenced under the [GNU GPL](https://www.gnu.org/licenses/). 

# How to cite fit3Dcspline
If you use fit3Dcspline to process your data, please, cite our [paper](https://www.nature.com/articles/nmeth.4661):
  * Yiming Li, Markus Mund, Philipp Hoess, Joran Deschamps, Ulf Matti, Bianca Nijmeijer, Vilma Jimenez Sabinina, Jan Ellenberg, Ingmar Schoen, Jonas Ries.  Real-time 3D single-molecule localization using experimental point spread functions. Nat. Methods 15, 367–369 (2018).
