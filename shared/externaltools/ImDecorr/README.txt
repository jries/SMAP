A detailled description of the method can be found in :+ 
Descloux, A., K. S. Grußmayer, and A. Radenovic. Parameter-free image 
resolution estimation based on decorrelation analysis.
Nature methods (2019): 1-7.

If this software has been usefull to your research, please consider citing

https://doi.org/10.1038/s41592-019-0515-7

-------------------
This README.txt file should be accompanied by

* The plugin
ImageDecorrelationAnalysis.jar in ijplugin folder

* A demonstration dataset
test_image.tif
This data set consists in a confocal and a STED image of COS7 microtubules labelled with alexa594 combined in a single stack

* A copy of the licence
LICENSE.txt

* A Matlab script
main_imageDecorr.m

* Matlab functions in a folder required by main_imageDecorr.m for image decorrelation analysis
funcs

-------------------
System requirements

ImageDecorrelationAnalysis_plugin.jar will require a running installation of ImageJ
Full instruction for the installation of ImageJ can be found at:
https://imagej.nih.gov/ij/download.html

The plugin requires ImageJ 1.48 or higher (tested for 1.48 to 1.52)
ImageJ version number can be found under: Help -> About ImageJ...

Matlab functions were developped using 64bits Matlab R2017b

-------------------
Installation guide

Copy the file ImageDecorrelationAnalysis_plugin.jar in your ImageJ plugins folder (.../ImageJ/plugins)
Restart ImageJ and the plugin will appear in the plugin list, under: Plugins -> ImageDecorrelationAnalysis

-------------------
Demo

Open the image of interest in ImageJ and make sure it is the active window. Before running the analysis, it is important to set the image pixel size and units (Image -> Properties… -> Unit of length, Pixel width and pixel height). If empty or not defined, the plugin will return the resolution in pixel.
Click on compute and once the analysis is done, results are added to the main results table.

The plugin supports multidimensional stacks of any bit depth. If a RGB image is supplied, it will be automatically converted to grayscale. The plugin also supports rectangular ROI, which allows to repeat the analysis on sub-regions and check the consistency of the estimate over the whole image.

The settings panel is composed of 4 optional parameters of the computation, specifying the range of normalized frequencies where the curve has to be computed (from Radius min to Radius max), as well as the number of points in between (Nr, typically 50 points). Finally, Ng (typically 10) specify the number of intermediate high-pass filtering used to find the resolution.

We also provide additional processing options in the form of 3 check-boxes:

*Do plot:
If checked, plot all the computed curves and results.

*Batch stack:
If checked, process all the frames, slices and channel of the current image. If do plot is also checked, all the curves for all the images will be plotted.

*Batch folder:
If checked, the user will be asked to select a folder containing images. All images are then opened, processed and the result table is automatically saved in the selected folder as a .csv file. Again, if “do plot” is checked, all the curves will be plotted.

-------------------
Instructions for use

1. Open ImageJ

2. Go to Plugins -> ImageDecorrelationAnalysis
A window will open with default settings already set

3. Open the demonstration data set
demo_COS7_a-tub_abberior_star635_confocal_STED.tif

4. Click on Compute

5. After about 15 seconds, the results of the analysis on the image will be added to the main results tabel and a figure showing the decorrelation curves will be generated and displayed

