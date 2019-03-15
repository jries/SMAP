# User Guide for SMAP
## 1	Overview
* Data vs. functionality (modules: basic modules and plugins) vs. GUI
* Define Module
* Difference to plugin

## 2	Data formats
### 2.1	Overview

SMAP has to process many different kinds of data. Examples include: a) raw camera images including meta data; b) the single-molecule localizations, which is a list of positions and additional properties such as number of molecules, size of the PSF etc.; c) reconstructed superresolution images and d) data from additional analysis such as statistics on the localization data, line profiles, cluster analysis,…; d) various parameters that are needed for the different processing steps (e.g. pixel size of reconstructed image, path to camera images,...). In the following paragraphs we describe how such data is represented in SMAP and how SMAP interacts with such data.

### 2.2	Camera images

The raw data of SMLM consists of a sequence of images containing the blinking single molecules. The pixel intensities are gray scale in Analog Digital Units (ADU). The acquisition parameters are usually saved as metadata and required to convert the pixel intensities from ADU to photons. We treat the metadata as parameters.

### 2.3	Localization Data

#### Fields

The single-molecule localization algorithm extracts the positions of the emitters from the images. Every single emitter, called ‘localization’ is then characterized by their x, y and z position and various other properties such as brightness, background, frame or localization precision. All these properties we call ‘fields’. You can picture the localization data as a huge table, where each row corresponds to a single localization and each column corresponds to one field. This localization data is the basis for further analysis such as reconstruction of superresolution images. Internally, every field is stored as a vector with a length corresponding to the number of localizations. 
There are two natural length units for localization data: the camera pixel size and nanometers. SMAP converts the units at the beginning to nanometers and the unit is usually apparent from the field name (‘xpix’ vs ‘xnm’).

#### Grouping

As the single-molecule switching is a stochastic process, many fluorophores emit not only during a single frame, but during multiple consecutive frames. As they belong the same event, often it is desirable to combine such localizations into a single localization. This is called grouping. However, some analysis and visualization algorithms work better on ungrouped data. As grouping has considerable computational complexity, in SMAP the grouping is not performed on the fly, but stored as a second localization data set. This allows rapid switching between grouped and ungrouped data.

#### Localization filters
Filters on the localization data select a sub-set of the localizations based on their fields. Filtering is an important step during data processing to a) remove bad localizations (e.g. low localization precision or photon numbers, high background, bad convergence of the fitter) and b) restrict the localizations to a certain spatial region. In practice, a variety of additional filters are used to select individual files, channels or frames or to perform sectioning in z by filtering the size of the fitted PSF and many more. As most analysis modules apply many different filters, SMAP pre-calculates the filters for each filtered field, re-calculates them when needed and combines them. This turned out to be faster than filtering on the fly.

### 2.4	Superresolution images

The most common visualization of single-molecule localizations is a rendered superresolution image. The pixel size of the superresolution image can be freely chosen. There are many different approaches to the rendering, with own advantages and disadvantages:

#### Scatter plot

Every localization is plotted as a dot.

#### Histogram rendering

Every pixel in the superresolution image has a value corresponding to the number of localizations falling into that pixel. For small pixel sizes at maximum one localization is found per pixel and histogram rendering visually approaches the scatter plot.

#### Gaussian rendering

Each localization is represented by a Gaussian intensity profile. In the superresolution image profiles of all localizations are added up. Usually, the integral of the Gaussian of each localization is one and the width of the Gaussian (sigma) is proportional to the localization precision. For fast rendering, often one sigma is chosen for all localizations. Then the superresolution image can be calculated from the histogram rendering by Gaussian filtering (blurring). SMAP can render each localization with a sigma proportional to its localization precision. A large proportionality factor leads to smoother details, but effectively reduces the resolution; a small factor avoids blurring of details but leads to more pointillist images. 

#### Diffraction limited rendering

By rendering each localization with a Gaussian corresponding to the size of the PSF and by setting the pixel size of the reconstructed image to the camera pixel size an image can be rendered which looks like a diffraction limited image. This can be useful for demonstration purposes, but the images can differ considerably from real diffraction limited images and usually show a higher contrast, since the background does not result in localizations.

#### Other rendering modes

Includes Quad-Tree Based Adaptive Histograms or Delaunay Triangulation Based Visualizations{Baddeley:2010hz}. In the end, it depends on the further analysis of the SR images and the preference of the researcher, which reconstruction method to use.
SMAP supports histogram and Gaussian rendering and rendered and original diffraction limited images and a modular architecture to easily incorporate additional reconstruction algorithms.
SMAP allows the overlay of many independently rendered images and external tiff files. Each of these images is referred to as a ‘layer’. Each layer has its own localization filters and specific parameters.

### 2.5	Further analysis

Rendering of a SR image is only one way of obtaining insights into a biological question. Further analysis can be based on the reconstructed SR image, but also be based on the localizations. The resulting data can be of various types and is usually very specific to the analysis algorithm.

### 2.6	Parameters

Most modules in SMAP depend on specific parameters. For instance, Gaussian reonstructor needs to know not only the single-molecule positions and localizations, but also the pixel size of the reconstructed image, its size in x and y, the look-up table to use etc. Many parameters belong to a specific module, but others are shared among modules. In SMAP, each module provides a description of a user interface that provides the local parameters. Global parameters are shared by all modules and can be set and read using a unique name. Fields of a user interface can be linked to a global parameter and thus linked among modules. User Interface

### 2.7	Overview

The graphical user interface (GUI) of SMAP provides easy access to all functionality to easily and intuitively analyze and explore the data. Upon startup it is in a compact form, but individual parts can be detached and positioned elsewhere on the screen (see customization). An overview of the main components can be found in Fig. xx.

### 2.8	Components

#### GUI for main functionality and plugins

Due to its extensibility SMAP consists of many modules. To easily find and use a specific module, the modules are arranged and ordered in the GUI in several levels. The main part of the GUI consists of a tab group, and each tab corresponds to a general class of functionalities (e.g. there is one tab for loading and saving data, one tab for localization of single molecules, one tab for rendering and several tabs for further data analysis).
In many of the tabs the individual modules are arranged again in a tab group or a list.
SMAP has 4 pre-defined tabs with custom functionality (File, Localize, Render, Siteexplorer) and 2 tabs for plugins (Process, Analyze). 

#### GUI for SR reconstruction settings

The dimensions of the reconstructed superresolution image are defined in SMAP by its pixel size, its center position (in nm) and its size (in nm). In this GUI global parameters for the reconstruction are set (e.g. pixel size). Additional parameters (e.g. size of reconstructed image and binning of pixels) can be accessed via the [Par] button. The [Reset] button resets the pixel size to fit the whole data set into the window. If not size is set in [Par], each pixel is rendered as a screen pixel and the size of the superresolution image is determined by the size of its window.
In addition a ROI can be associated to the image. SMAP allows for point [+], square [ [] ], elliptical [O], poligonial [<>], free [{}] or line [ | ] ROIs. For line ROIs the width of the roi is specifited in the field |<->|. The show checkbox toggles re-drawing of the ROI.

#### Localization filter GUI

This part of the GUI performs filtering of the localization data. Every layer has its own instance of this GUI. It consists of two windows. The upper window shows a table in which every row corresponds to one field of the localization data. The minimum, mean and maximum for each field are calculated. In addition, the minimum and maximum can be set to filter on the corresponding field. For filtering, the filter checkbox needs to be checked.
Selection of a specific field opens its histogram view in the window below. Here the values of the selected field are displayed as a histogram. The black histogram corresponds to all values, the green histogram to the selected values and the red histogram shows the values of the filtered localizations shown in the specific layer. 
Functionality of the histogram gui: Sliders and edit fields for the minimum and maximum value, filter on/off, auto-update (if selected, the SR image is reconstructed each time a value is changed), range fix (fixes the range to the specified parameter when moving the sliders).
A subset of the fields can be also selected in the Render/Layer tab.
This GUI might be hidden by the overview image. To make it visible, use the [OV -> filter] button in the settings GUI.

#### Overview image

The overview image shows a reconstruction of the whole data set. The size of the SR reconstruction is denoted by a box. [Update] re-renders the overview image. Clicking on the overview image centers the SR image on this point.

#### Menu

The menu provides quick access to functionality. It has a general settings item to define parameters which are loaded on startup. The File menu provides access to load and save functionality (also found in the main GUI/File tab), the plugin menu provides access to all plugins available in SMAP, organized by their directory structure. Additional menu items can be user defined as described in customization.

#### Other windows

SMAP opens additional windows such as the superresolution reconstruction, additional GUI components of plugins or the results of plugins.

## 3	Plugins
### 3.1	Modules in SMAP

### 3.2	Workflow plugins

### 3.3	Localization data plugins

## 4	Single-molecule localization

### 5	Superresolution rendering

### 5.1	Layers

## 6	Customization (move to end? After WF etc)

SMAP can be customized to provide fast access to only the used functionality without cluttering the GUI with unused functionality. Most customization options can be accessed with a context menu (right-click).

### 6.1	Detaching windows

The reconstruction settings GUI, the filter GUI, the overview image, any tab of the main GUI and any plugin can be detached and opened in a separate window. This makes them easily accessible on a big screen. Note that detaching is currently not possible. Just open SMAP again to have the compact GUI.

### 6.2	Main GUI items

The sub-tabgroups of the main GUI offer customization options. In the analyze and process tabgroups, new tabs can be created and tabs can be deleted. The new tabs can be populated with plugins (context menu of list of plugins). This allows full customization of the access to the plugins. Changes are automatically saved and loaded on restart and can be reset in the plugin list context menu and the menu. The guimodule settings are saved in ‘settings/temp/guimodules.txt’. Deleting this file populates the GUI with all plugins in the processor and analyzer directories. Copy a specific configuration (guimodules.txt) to ‘settings/guimodulesdefault.txt’ to use this as a default. 

### 6.3	Localization GUI

### 6.4	Menu customization


## 7	Specific plugins and tasks
