---
title: 'SMAP Step-by-Step Guide'
---

Installation
============

Requirements
------------

1.  Matlab newer than 2014a, performance boost with 2015b. Some
    functions require 2016a or newer. Toolboxes: Optimization, Image
    processing, Curve fitting, Statistics and Machine Learning.
    Optional: Parallel Computing. A stand-alone version will be
    released, but will be limited in extensibility.

2.  Mac or Windows

3.  For GPU fitting: Windows, NVIDIA graphics card. CUDA driver
    (recommended: version 7.5).

Installation
------------

1.  Clone git repository:

    a.  Use Terminal (MacOS) or Cmd (Win). Use cd to navigate to the
        target directory. (e.g. cd git)

    b.  Type: `git clone` <https://github.com/jries/SMAP-light> and type
        in your username and password for your git account.

    c.  Install the 3D fitter by typing

> `git clone` https://github.com/jries/fit3Dcspline-light

2.  Install Micromanager 1.4.22 from https://micro-manager.org

3.  If needed install www.openmicroscopy.org/bio-formats/downloads

4.  In Matlab: run SMAP.m, if questioned, change folder.

5.  In the Menu select SMAP/Preferences... Switch to the Directories tab
    and select the directories of Micro-Manager and of the
    bioformats\_package.jar. Press **Save and exit**.

Micromanager users
------------------

1.  Save your images as single or multi-image Tiff stacks, turn on "save
    metadata" in preferences.

Camera Settings
---------------

> Add your camera and camera settings to SMAP

1.  Acquire a data set with the right camera settings.

2.  Open Menu: SMAP/Camera Manager

3.  Use **Load images** to load the data set. If the camera has not been
    recognized you will be prompted to add a new camera. Otherwise you
    can add and remove cameras with a right click on the camera list.

4.  The default camera can be used if no camera is recognized.

5.  The camera manager allows extracting acquisition parameters from the
    metadata. SMAP uses the following parameters:

    a.  EMon (logial, if gain is on)

    b.  Cam\_pixelsize\_um (pixel size in the object plane in
        micrometers)

    c.  Conversion (in ADU/e-, used to calculate photons from the camera
        units, by default it does not include the em-gain)

    d.  Emgain (em gain value)

    e.  Offset (camera image offset in ADU)

    f.  Roi (ROI on the camera chip in pixels)

    g.  Exposure, timediff, comment (for information only)

    h.  numberOfFrames (number of frames)

    i.  Width, Height of the image

6.  For each parameter that SMAP uses define the mode:

    j.  fix: uses the value in the column 'fixvalue'.

    k.  metadata: uses the field defined in 'metafield' together with
        the parser defined in the 'conversion' column. Here *X* is
        substituted by the metadata corresponding to the metafield.
        Choose the metafield by clicking on it.

    l.  State dependent: This allows you to define camera parameters in
        dependence on camera settings. This is useful e.g. for the
        conversion factor. The state is defined by the parameters in the
        list *state defining parameters*. Select metadatafields, those
        define together with the concrete values the camera state.
        Define the values of the SMAP parameters in the list to the
        right.

Single-molecule localization
============================

Selection of a localization Workflow
------------------------------------

1.  A workflow for single-molecule localization can be selected in the
    **Localize** tab with the **Change** button close to the bottom
    right of the window. For most cases 'fit\_fastsimple' is a good
    choice.

2.  You can find a description of the current fitting workflow and a
    graphical representation of the modules by pressing **Info** next to
    the **Change** button.

Basic fitting
-------------

1.  In the **Localize** Tab, **load images**: Select one image inside a
    directory containing all the tifs, or a tiff-stack. Alternatively,
    you can select any OME-readable file, but import of metadata then is
    limited.

2.  If sufficient metadata could be extracted and if the camera was
    added in the Camera Manager, the acquisition parameters are
    automatically set. Otherwise you can either **load metadata** from a
    previous experiment or manually **set Cam Parameters**.

3.  You can specify a frame range which to fit. Check **Online
    analysis** if you want to fit during the acquisition (then the
    maximum frame is ignored, and SMAP waits for new images).

4.  In the **Peak Finder** Tab you can set the parameters for the
    initial guessing of single-molecule positions. Usually this is done
    on a background-corrected image or on an image convoluted with a
    *difference of Gaussians (DoG)*. Use ToolTips (hover mouse over
    control) to get information about specific parameters.

5.  By pressing **Preview** after selecting a frame with the slider next
    to it an image will open which will show the positions of the found
    localizations. Use this to optimize peak finding parameters. The
    preview mode determines which images are shown here.

6.  You can restrict the fitting to pre-defined regions in the image.
    First select either an elliptical or rectangular ROI in the popup
    menu. With **Fit in ROI** you can select a region in which to fit,
    with **Remove ROI** you can select a region which to exclude from
    the fit. You can define multiple ROIs. **Clear ROI** to fit
    everything.

7.  In the Fitter Tab you set

    a.  The size of the ROI in which the fitting of single molecules is
        performed

    b.  The fitter module and its parameters

    c.  If to fit on the background corrected data (not recommended)

8.  In the **Localizations** Tab you can switch on the rendering during
    the fitting and set the update time. In addition, you can do simple
    pre-filtering of the data. To set these parameters you can again use
    **Preview**. Found candidate positions are marked with a box, found
    localizations which pass the filters are marked by a circle.

9.  Now you can fit the whole data by pressing **Localize**.

10. The fitted localizations are automatically saved in the base
    directory of the images with an extension '\_sml'.

Batch processing
----------------

1.  You can save your acquisition fit settings (previous paragraph) with
    the **Batch** button.

2.  In the **Input Image** Tab open the **Batch processor**.

3.  The batch file you just saved is already set as the main batch file.
    But you can replace it by another with **load main batchfile.** If
    **use for all** is checked, this will be used for all the fits,
    otherwise only for the datasets which are not imported to the batch
    processor with a batch file.

4.  With **add** you can add a) further batch files, b) one single image
    from a stack or c) a multiple tiff stack. These appear in the list
    on the left.

5.  You can multiple directories with **add directories**. These
    directories contain a) tiff images, b) further directories with Tiff
    images inside (here use the filter string below to specify which
    directories to load, and the **\>\#images** to set a lower limit for
    the number of images.

6.  You can remove items from the list, but don't empty it. With Batch
    process the fitting starts.

7.  If you **add online directory** and start the batch processor, it
    checks for new fittable directories in this directory and
    automatically fits them (used e.g. for automated microscopy).

Rendering
=========

Load localizations
------------------

1.  In the **File** tab press **Load** and select a file containing
    localizations ('\_sml.mat', but also '.csv').

2.  The localizations are automatically grouped (i.e. localizations in
    adjacent frames stemming from one and the same fluorophore are
    combined into one localization), using the parameters below (maximum
    allowed displacement, maximum time in frames the molecule can be
    dark). Note that also the ungrouped localizations are always
    available. If you change the parameters, press **Group** to regroup.

3.  After loading, the **Render** tab is opened and an overview image is
    displayed. By clicking in the overview image or pressing
    **Reconstruct**, the superresolution image is calculated.

4.  **Load** clears current data before loading. **Add** adds a file to
    the current localizations without clearing already loaded
    localizations.

5.  After loading localizations, you can add single tiff images
    (diffraction limited markers) and associate them to a localization
    file.

Modify the size and location of the image 
------------------------------------------

1.  Set the pixelsize in the Format GUI (or use the mouse wheel to zoom
    in and out). Use pre-defined pixel sizes.

2.  You can change the size of the image window.

3.  You can move around the superresolution image by clicking on it,
    then the clicked point will be centered.

4.  Right-click resets the view to display all localizations.

5.  In **par** you can specify a size of the image different from the
    screen resolution and binning of pixels.

6.  Here you can also quickly turn individual layers on and off.

Set the appearance of the image
-------------------------------

1.  The parameters in the **Layer** tab determine how the image is
    rendered and which localizations are rendered.

2.  You can define multiple layers (**+**), which are overlaid in the
    superresolution image.

3.  The checkbox in the upper left corner determines if the layer is
    displayed.

4.  Select the file and the channel(s) to display.

5.  Select the renderer: 'Gauss', 'histogram' 'constant Gauss' or
    'diffraction limited' reconstruction or 'Other' (external renderer).
    In case you have attached a Tiff image to the data, you can choose
    here to display it.

6.  You can determine the integrated intensity of single localizations.
    It can be equal to 1 (normal), equal to the number of photons or
    equal to the number of blinks (only makes sense if grouped data is
    displayed).

7.  Select Color-coding: 'Normal' uses the value of the reconstructed
    image for coloring, but you can also color the image according to
    the z-coordinates or any other field (property) of the
    localizations.

8.  Select the corresponding look-up table (LUT).

9.  The values of c-range determine the range of the parameters used for
    coloring that are mapped onto the entire LUT. Use **remove out** to
    remove localizations outside the LUT, otherwise they will be
    displayed with the minimum or maximum color of the LUT.

10. You can select with the **group** checkbox if to display grouped or
    ungrouped localizations.

11. The contrast button brings up a histogram to select which image
    intensity is mapped to the maximum intensity (higher intensities are
    saturated). You can toggle to use **absolute** intensities, or the
    fraction of pixels to be not saturated (**quantile**). The quantile
    parameter can be between 0 and 1 (typically: 0.995) or a negative
    number Q (typically -3.5). Then the fraction 10\^Q is not saturated.

12. With the remaining fields you can determine minimum and maximum
    values for filtering (see paragraph below).

13. The **par** button opens a dialog to set additional parameters:
    minimum size of the Gaussian image for reconstruction, the size of
    the reconstructed Gaussian in units of the localization precision.

14. Shift x,y shifts the image in the associated layer, this can be used
    to correct for shifts between images.

15. By pressing **save** you can save the current settings of the Layer.
    Pressing **default** loads those settings. **-\>all** copies the
    settings of the current layer to all other layers.

Filtering of localizations
--------------------------

1.  You can toggle between the overview image and the filter GUI by
    pressing **OV-filter**.

2.  The upper table lists all properties (fields) of the single molecule
    localizations together with their minimum, mean and maximum value.
    You can set minimum and maximum values. Importantly, you can select
    if to filter these on these fields with the checkbox.

3.  Below, you see a histogram representation of a specific field. You
    can select the field by either pressing on a row in the table or on
    a button in a **Layer** corresponding to the specific field (locp,
    frame, PSF, locprec z, z)

4.  You can switch the filter on and off and change the range with the
    sliders. If **Auto update** is checked, the image is directly
    rendered on the fly. If you check **range fix** and move the
    sliders, the difference between minimum and maximum slider is fixed
    to the value below **range fix**.

ROIs
----

1.  You can define a region of interest with the buttons in ROIs (Format
    GUI). These ROIs are used by various plugins. For Line-ROIs you can
    specify the thickness of the ROI. You can toggle redrawing of ROIs.

Saving
------

1.  In the **File** tab you can select what to save (localizations,
    settings, Tiff images) and press **save**.

2.  When saving localizations (as '\_sml.mat' or '.csv') you can check
    **only visible** to save only the localizations currently displayed.

3.  When saving Tiff images, the reconstruction parameters are saved as
    comments. When opened in Fiji, the size shown in the title is in
    nanometers (although Fiji calls them pixels).

ROI Manager
===========

The ROI manager allows for simple automated, semi-automated or manual
selection of ROIs that can be then annotated, sorted and run through an
evaluation pipeline. The results of this evaluation can then be
statistically analyzed.

Manually generating a list of ROIs
----------------------------------

1.  Tab **ROIs/Settings** click **show ROI manager**

2.  You find panels that show the superresolution image of a whole file,
    a part of the file (called cell) or a ROI are, as well as lists to
    select stored ROIs and cells. The file list is linked to the file
    list of the **File** Tab and cannot be edited here.

3.  You can define the pixelsize for reconstruction and the FoV for the
    cells and sites (regions around the ROIs), as well as the ROI size
    itself.

4.  Check **rotate** to rotate sites and **draw boxes** to show the
    positions of the selected cells or sites.

5.  By clicking on an item in a list in the **Roi manager** you can
    select and draw it.

6.  For fast scrolling through sites, the reconstructions are saved.
    Therefore, if you change any parameters (e.g. size of the FoV, or
    render parameters in the **Layers**) you need to **redraw**. You can
    **redraw all** in the **ROI/Settings** Tab.

7.  By left-clicking in the File image you can select a cell. Add it
    with the **+** button above the cell image to the list. You can move
    a cell by right-clicking in the cell image. The cell will be
    centered on that spot.

8.  In the same way you can left-click in the cell image to define a
    site and add it with the **+** button, and move it by right-clicking
    in it.

9.  You can rotate a ROI by pressing **Angle** and drawing a line. The
    ROI is rotated so that the line is horizontal.

Annotate ROis manually
----------------------

1.  In the tab **ROIs/Annotation** ROIs can be manually annotated. There
    are four lists to choose from. The items are defined in text files
    (look at 'settings/parlistdefault.txt and modify accordingly). You
    can load your settings file with **\<-load**.

2.  Use left-arrow and right-arrow keys to go to the previous and next
    site, respectively

3.  Annotate site by clicking on the lists. Keyboard shortcuts are:
    up-arrow and down-arrow to choose list entries and shift +
    left-arrow/right-arrow to go to the previous/next list.

4.  You can draw additional lines (two buttons on the right) and
    annotate size and angle.

5.  Add an additional comment if needed.

Sorting of ROIs
---------------

1.  In the tab **ROIs/Sort** you can sort the ROIs according to up to
    four criteria.

2.  Define if you want to sort ascending or descending

3.  Select the parameter that is used for sorting:

    a.  Hierarchy: File, Cell, Site

    b.  Statistics: Number of photons, PSF size etc...

    c.  List: any of the lists

    d.  Annotation: any of the lines (length)

    e.  Evaluation: the results of any evaluator. You can choose it with
        the **select** button. Use the list to navigate through all
        results.

    f.  Other: you can select any parameter saved with a ROI for
        sorting.

4.  Sort with the **Sort** button.

Evaluation
----------

1.  In the tab **ROIs/Evaluation** you can select several evaluation
    processors, which evaluate each site and return results that are
    then saved with the ROIs.

2.  Select processors with **add module** and **remove**.

3.  If a module is checked, it is used for evaluation.

4.  Clicking on a module in the list opens its GUI on the right. You can
    adjust parameters here and run the evaluation with **preview** or by
    redrawing a ROI.

5.  Re-evaluate all ROIs with the same settings with **redraw all** in
    the **ROIs/Settings** tab.

Analyze evaluation results
--------------------------

1.  In **ROIs/Analyze** you can find plugins to analyze results.

Automatic segmentation
----------------------

1.  In **ROIs/Segment** you can find plugins to automatically segment
    files and store the result as ROIs.

Selected plugins
================

1.  Plugins are found in the **Plugins** menu. A selection of regularly
    used plugins can be found in the **Analyze** and the **Process**
    tabs in subtabs (configurable).

2.  Select a plugin, edit the parameters and press Process.

3.  With **showresults** you can toggle the window with the output of
    the module on and off.

4.  **Info** displays a description text of the module in the results
    window.

Drift correction 
-----------------

*Process/drift/driftcorrection*

Drift correction based on the localizations, but works also very well in
case fiducial markers are present (in that case render the image
ungrouped).

1.  Select parameters to render a large part of the image. Only the FoV
    of the superresolution image is used for drift correction.

2.  Choose the number of time points to perform the drift correction on
    (typically 7-25, this algorithm rather corrects for drifts than for
    fast jumps or oscillations). The other parameters usually need not
    be optimized (use Tool Tips to understand what they mean).

3.  Use **Reference is last frame** to drift correct the first of two
    consecutive measurements.

4.  Press **Run**. With show results you can display the results of the
    procedure.

5.  The drift-corrected localizations are automatically saved as
    '\_driftc\_sml.mat' files.

Localization statistics
-----------------------

*Analyze/measure/statistics*

Get single-molecule statistics

1.  If **use Roi** is checked, only the localizations in the current ROI
    /FoV are evaluated.

2.  If **use layers/filters** is checked, each layer is evaluated
    individually; otherwise statistics for grouped and ungrouped data
    are shown.

3.  With **plot overview** you can have all results in one figure (e.g.
    for saving) rather than in individual tabs.

FRC resolution
--------------

*Analyze/measure/FRC resolution*

Calculates the FRC resolution according to: R. P. J. Nieuwenhuizen, K.
A. Lidke, M. Bates, D. L. Puig, D. Grünwald, S. Stallinga, and B.
Rieger, "Measuring image resolution in optical nanoscopy," Nat Methods,
vol. 10, no. 6, pp. 557--562, Apr. 2013.

3D viewer
---------

*Analyze/sr3D/Viwer3D*

1.  Define a linear ROI in the superresolution image.

2.  Press **Process**. A side-view reconstruction is opened.

3.  You can **set pixelsize** manually; otherwise the pixel size of the
    current reconstruction is used.

4.  You can use the controls to translate, rotate or zoom. '0' resets
    the view.

5.  When the sideview window is selected and on top, you can use key
    shortcuts to translate, rotate (command / strg) or zoom (alt, this
    changes the size of the ROI). The direction is defined by the arrow
    keys. The direction perpendicular to the screen can be accessed with
    the '.' and ',' keys.

6.  Pressing 'shift' results in a smaller movement.

7.  You can also manually move the ROI in the superresolution image, the
    3D reconstruction is updated on-the-fly.

8.  With **use transparency** localizations closer to you partially
    block localizations in the background for a better 3D look.

Calibrate Astigmatic or any other 3D PSF
----------------------------------------

*Analyze/sr3D/calibrate3DsplinePSF*

This plugin generates a cspline model of your experimental PSF (see Li
*et al,* 2018. Real-time 3D single-molecule localization using
experimental point spread functions. Nature Methods). For this you need
to acquire several z-stacks of beads immobilized on a coverslip. A range
from -1 µm to 1 µm with respect to the focal plane and a distance
between the planes of 10-50 nm works well. Please consult the
*User\_guide\_Ries* in the fit3Dcspline-light directory for further
information.

1.  Press **run** to open the GUI to calibrate the 3D PSF.

2.  **Select camera files** to open a file dialog box to select several
    files at the same time

    a.  With **add** you can add a single file or several files in the
        same directory

    b.  With **add dir** you can add several directories. SMAP will
        automatically try to find image files in those directories and
        add them.

    c.  Press **Done**

3.  The output file is set automatically, but you can change it manually
    with **Select output file**.

4.  Select the 3D modality: *arbitrary* is the most common choice.

5.  If your PSF does not show strong variations along *z* or a symmetry
    with respect to the focus (e.g. because it is an unmodified 2D PSF)
    check **2D.**

6.  Select the distance between the frames you used for acquiring the
    z-stacks.

7.  Leave other parameters as they are, but if you have dense beads
    decrease **Minimum distance** and **ROI size**.

8.  **Calculate bead calibration** calculates and saves the PSF model

You can use this model for fitting as described above using e.g. the
*fit\_fastsimple* workflow using the following modifications:

9.  In the *Fitter* tab select **Spline** as the fitting model.

10. Load the PSF model you previously generated with **Load 3D cal**.

11. If you have a symmetric PSF you can use different start parameters
    that you can insert in **z start (nm)**. Default is **0**.

Mathematics
-----------

MathParser

Multi-color and registration
============================

Find a transformation to register two data sets
-----------------------------------------------

*Process/Register/Register Localizations*

1.  Select the target (usually bottom or right), and if to mirror the
    target half-image.

2.  Select a transformation type (try projective, if that is not
    sufficient use polynomial with a parameter 3).

3.  Under Parameters you can choose additional parameters

    a.  Pixel size for correlation. Around size of the localization
        precision. If the correlation image is dotty and the wrong
        maximum is found, increase this size.

    b.  Max shift for correlation: reduce, if wrong maximum is found.
        Increase, if true maximum is outside.

    c.  Max locs for matching: eg. 100 000. Numer of localizations used
        to determine transformation. Precision increases with this, so
        does computation time.

    d.  Max shift matching: distance that corresponding localizations
        can be apart (after shift is applied). 250-500 nm typically. If
        this value is too large, random localizations are matched, this
        can introduce systematic error.

4.  Press **Run** and judge results:

    e.  The *shiftcorr* should show a clear maximum, the square should
        be on that maximum. If you see many dots around this maxiumum,
        increase the pixel size.

    f.  The *scatter* should show a clear maximum in the center. *dx*
        and *dy* should be 20-80 nm. The number of anchor points should
        be at least a few %. *hist* is just a profile through the
        scatter image.

    g.  If the transformation is good, save it with **save T**. Note
        that in other plugins the default localization file is
        initialized with this file, making it optional to load a
        transformation file.

5.  For difficult data you can also first find an approximate
    transformation (e.g. projective) as described before (or load one,
    with **load T**). Then check **use initial T**. This transformation
    is then applied before finding anchor points.

6.  If **use layers** is checked, the plugin does not use all
    localizations, but only those displayed in T: and R: (e.g. for two
    synchronized cameras, then use center for target).

Register Images
---------------

Get intensities from camera images
----------------------------------

*Process/Assign2C/2C intensities from images*

This plugin uses a transformation (determined e.g. with *register
localizations*) to find for every localization the position in the other
channel and then determines the intensity in both channels.

1.  Load a transformation

2.  Per default, this plugin does median filtering. Select the spatial
    and temporal spacing for this (dx, dt).

3.  Select one or several plugins which determine the intensity:

    a.  Roi2int\_sum: uses a ROI (set size) to determine intensity, and
        a larger ROI for the background.

    b.  Roi2int\_fit: Uses a Gaussian fit to determine intensity and
        background. The position is fixed to the fitted position. You
        can use the fitted PSF size or fix it. If **fit on BG** is
        checked, the background is subtracted prior to fitting and the
        fit is performed with background set to zero. Otherwise the
        background is a fitting parameter.

4.  Press **Run** and when asked select the original raw camera images.
    The results are automatically saved with the \_dc in the file name.

Determine channel from intensities
----------------------------------

*Process/Assign2C/ Intensity2Channel*

This plugin assigns a channel value to the localizations based on two
fields of the localization data (usually intensity in camera channel 1
vs camera channel 2).

1.  Select the two fields in **value 1** and **value 2**.

2.  Press **Run** to display a 2D histogram (normal rendering and
    logarithmic color rendering which shows better any background).

3.  Adjust all parameters to obtain an optimal separation of the two
    clouds belonging to both dyes.

4.  Dye 1 and dye 2 are assigned channels 1 and 2, respectively.
    Localizations thich are excluded (blue in the *log split* tab) are
    assigned channel 0.

Apply Transformation
--------------------

*Process/Register/Apply Transformation.*

This plugin applies a saved transformation to localizations or images.
You can use this to move all localizations from the second channel into
the first channel.

1.  Load a transformation with **load T**.

2.  Select a dataset. If **transform all files** is selected, all loaded
    files are transformed with the same transformation.

3.  Select what to transform: reference (to target), target (to
    reference) or all (using either reference to target or target to
    reference transformation).

4.  Select if to transform only localizations or tiffs or both.

5.  If **set channel** is selected, the channel field is overwritten
    depending on the localizations being reference or target
    localizations.

6.  **Run**.

Combine Channels
----------------

*Process/Register/Combine Channels.*

Similar to *Apply Transformation.* In addition, fluorophores in both
channels which can be associated to a single fluorophore (same frame,
close proximity after transformation) are averaged and presented as a
single localization. This is useful for e.g. ratiometric dual-color
imaging or bi-plane 3D imaging to avoid duplication of localizations.

Ratiometric Dual-Color Imaging
------------------------------

1.  Fit your data

    a.  2D or 3D fit.

    b.  Do not yet apply drift correction. This you can do later after
        channel assignment.

2.  Find transformation for both channels with
    **Process/Register/Register Localizations.**

3.  Determine intensities of localizations in both channels with
    **Process/Assign2C/2C intensities from images.**

4.  Assign channel from relative intensities with
    **Assign2C/Intensity2Channel.**

5.  You can now render both channels individually with in two layers.

6.  Optional: transform target localizations onto reference with
    **Process/Register/Apply Transformation** or
    **Process/Register/Combine Channels**.
