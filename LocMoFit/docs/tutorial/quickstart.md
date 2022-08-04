# Quick start

In this tutorial you will learn how to perform your first model fitting with *LocMoFit*, in SMAP with a GUI. We will fit a 2D ring model to the top-view projection of NPCs. You will need:
* Software: **SMAP** installed. Further information can be found on our [GitHub](https://github.com/jries/SMAP/) site.
* Localization data: _U2OS_Nup96_BG-AF647_demo_sml.mat_

The mat file can be downloaded [here](https://www.embl.de/download/ries/LocMoFit/).

## Preparation
1. Start **SMAP** ({doc}`how to? <../howto/SMAP.runSMAP>`).
2. Load the localization data ({doc}`how to?<../howto/SMAP.loadData>`) _U2OS_Nup96_BG-AF647_demo_sml.mat_. This file contains segmented nuclear pore complexes that you will be analyzing in the following steps.
:::{note}
You can check out all the pre-processing steps from _fitting raw data (raw camera frames)_ to _segmentation of NPCs_ in [**SMAP_manual_NPC.pdf**](https://www.embl.de/download/ries/Documentation/SMAP_manual_NPC.pdf).
:::

## Loading LocMoFit
:::{include} ../howto/SMAP.loadLocMoFit.md
---
start-line: 1
---
:::
:::{note}
The button **?** at the up-right corner of the LocMoFit GUI provides you details of fields and buttons in the current tab/sub-tab.
   ![Help button](../images/helpButton.PNG)
:::

## Setup
We will be using the in-built model _{class}`ring2D<models.ring2D>`_ here ({doc}`more about model types<../basics/geometricModel>`). The following steps show you how to load the model into LocMoFit:


1. In the right panel, go to **[M1]** -> **[Model]**, click the drop-down menu (_arc2D_ is shown by default), and then select _ring2D_. Here _M1_ stands for _model 1_.
   ![Model tab](../images/modelTab_default.PNG)

2. Click **load model**. Now the model is loaded and LocMoFit is ready to fit.
	:::{Note}
	Click the button **i** next to the drop-down menu opens the webpage detailing the selected model.
	![Model info](../images/button_modelInfo.PNG)
	:::
## Fitting
Next, we will execute the fitting. This is done by clicking the site in the list of sites in the _ROI manager_ window:
1. Go to **[ROI]** -> **[Settings]**, click **show ROI manager**. This opens the **ROIManager** in a new window.
	![ROIManager](../images/ROIManager.PNG)
	:::{Note}
	**ROIManager** allows you to manage ROIs in different cells and files. Check section 8.2 (_Manually generating a list of ROIs_) in [**SMAP_manual_NPC.pdf**](https://www.embl.de/download/ries/Documentation/SMAP_UserGuide.pdf) for more information.
	:::
2. Go back to **[Evaluate]**, check **evaluate on** and **display** in the left panel.
	:::{Note}
	* The loaded modules will be evaluated only when **evaluate on** is checked.
	* Result windows of the loaded modules will be displayed only when **display** is checked.
	![Evaluate panel](../images/evaluatePanel.PNG)
	:::
3. In the _ROI Manager_ window, click on one site in the list of sites and wait for a few seconds. You should see a new window _LocMoFitGUI_, in which the localizations are plotted on the fitted model.
	![Viewer](../images/viewer_quickstart.PNG)
	:::{Note}
	The **LocMoFitGUI** window allows you to inspect the result right after fitting. The window may look different depending on the model type.
	:::
Now you have your first fit done! Congratulations! You can further explore a few sites to get familiar with the interface.

## Next tutorial
You are in the introductory series. The next tutorial is {doc}`Composite model<../tutorial/compositeModel>`.
