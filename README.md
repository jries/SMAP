# SMAP
Superresolution microscopy analysis platform

COPYRIGHT:     Jonas Ries, 2020
LICENSE:       GPLv3
AUTHOR:        Jonas Ries, EMBL Heidelberg, ries@embl.de 27.03.2020
               www.rieslab.de, www.github.com/jries/SMAP

Documentation and installation instructions at: 
    /Documents/Manual/SMAPStep-by-StepGuide.md
Instructions for using SMAP to analyze nuclear pore complex standard samples: 
    /Documentation/Manual/SMAP_manual_NPC.pdf

Requirements
------------

1.  Matlab 2016a or newer. Toolboxes: Optimization, Image
    processing, Curve fitting, Statistics and Machine Learning.
    A stand-alone version will be released, but will be limited in extensibility.

2.  Mac or Windows

3.  For GPU fitting: Windows, NVIDIA graphics card. CUDA driver
    (recommended: version 7.5).

Installation
------------

1.  Clone git repository:

    a.  Use Terminal (MacOS) or Cmd (Win). Use cd to navigate to the
        target directory. (e.g. cd git)

    b.  Type: `git clone` <https://github.com/jries/SMAP> and type
        in your username and password for your git account.

    c.  Install the 3D fitter by typing

> `git clone` https://github.com/jries/fit3Dcspline

2.  Install Micromanager 1.4.22 from https://micro-manager.org

3.  If needed install www.openmicroscopy.org/bio-formats/downloads

4.  In Matlab: run SMAP.m, if questioned, change folder.

5.  In the Menu select SMAP/Preferences... Switch to the Directories tab
    and select the directories of Micro-Manager and of the
    bioformats\_package.jar. Press **Save and exit**.

Typical install times: 15 minutes.