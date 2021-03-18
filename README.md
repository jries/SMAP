# SMAP
Superresolution microscopy analysis platform

COPYRIGHT:      Jonas Ries, 2020
LICENSE:        GPLv3
AUTHOR:         Jonas Ries, EMBL Heidelberg, ries@embl.de 27.03.2020
                www.rieslab.de, www.github.com/jries/SMAP
PLEASE CITE AS: Ries, J. SMAP: a modular super-resolution microscopy analysis
                platform for SMLM data. Nat Methods (2020).
                https://doi.org/10.1038/s41592-020-0938-1
                Please also cite the references for the plugins you use
                (as mentioned in the plugin info).

Documentation and installation instructions at:
    /Documents/Manual/SMAPStep-by-StepGuide.md
Instructions for using SMAP to analyze nuclear pore complex standard samples:
    /Documentation/Manual/SMAP_manual_NPC.pdf

Please cite as: Ries, J. SMAP: a modular super-resolution microscopy analysis platform for SMLM data. Nature Methods (2020) doi:10.1038/s41592-020-0938-1.


Requirements
------------

1.  MATLAB 2019a and newer. Toolboxes: Optimization, Image processing,
    Curve fitting, Statistics and Machine Learning.
    A fully functional stand-alone version that does not require a MATLAB
    license but is limited in extendability can be downloaded from https://www.embl.de/download/ries/SMAPCompiled/ (Installation notes see below).

2.  Mac or Windows

3.  For GPU fitting: Windows, NVIDIA graphics card. CUDA driver
    (recommended: version 7.5). All fitters also come with a CPU version
    that is used when these specifications are not met.

Installation
------------

1.  Clone git repository:

    a.  Use Terminal (MacOS) or Cmd (Win). Use cd to navigate to the
        target directory. (e.g. cd git)

    b.  Type: `git clone https://github.com/jries/SMAP` and type
        in your username and password for your git account.

    c.  Install the 3D fitter by typing
        `git clone https://github.com/jries/fit3Dcspline`

2.  Install Micromanager 1.4.22 from https://micro-manager.org

3.  If needed install www.openmicroscopy.org/bio-formats/downloads

4.  In MATLAB: run SMAP.m, if questioned, change folder.

5.  In the Menu select SMAP/Preferences... Switch to the Directories tab
    and select the directories of Micro-Manager and of the
    bioformats\_package.jar. Press **Save and exit**.

Typical install times: 15 minutes.

Installation of stand-alone version
------------

1.  Download the respective version from https://www.embl.de/download/ries/SMAPCompiled/ corresponding to your operating system.
2.  Follow the installation notes in Installation_notes_SMAP_compiled.rtf which can be found under the previous link.
