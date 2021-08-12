Welcome
=======

You can find all things related to LocMoFit here.

Overview
--------

LocMoFit is a model fitting framework for SMLM data.

Menu info.
----------

If you are new to LocMoFit, you can find out how to get started via __:doc:`Getting statrted`__ in the menu on the left. We recommend you to install LocMoFit with SMAP, a modular super-resolution microscopy analysis platform for SMLM data we developed.

.. important::
   We recommend you to start with our introductory series by following the tutorial in :doc:`Quick start`.

__:doc:`Structure`__ gives an overview of what LocMoFit has and the relations between different classes.

Pages under __BASICS__ contain the necessary information for working with LocMoFit.

Pages under __TUTORIAL__ get you familiar with LocMoFit through hands-on tutorials.

Pages under __PROGRAMMING NOTE__ contain information for working with LocMoFit through coding.

You can find out other resources LocMoFit provides in pages under __RESOURCE__ and API documentations in pages under __API__.


.. toctree::
   :maxdepth: 2  
   :hidden:
   
   start
   structure
   
.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Basics
   
   ./basics/quickstart.md
   ./basics/geometricModel.md
   SMAP.md

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Tutorial
   
   Composite model <./tutorial/compositeModel.md>
   ./tutorial/chainSteps.md
   ./tutorial/fourColorNPC.md
   
.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Programming note
   
   ./programming/buildModel.md
   
.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: API documentation
   
   modLibrary
   references
