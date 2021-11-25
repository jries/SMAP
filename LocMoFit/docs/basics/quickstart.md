# Quick start

In this tutorial you will learn how to perform your first model fitting with *LocMoFit*, in SMAP with a GUI. We will fit a 2D ring model to the top-view projection of NPCs. You will need:
* **SMAP** installed. Further information can be found on our [GitHub](https://github.com/jries/SMAP/) site.
* _U2OS_Nup96_BG-AF647_demo_sml.mat_
* _ring2d.png_

The two files above can be downloaded [here](https://www.embl.de/download/ries/LocMoFit/).

## Preparation
1. Start **SMAP**.
2. Load the dataset _U2OS_Nup96_BG-AF647_demo_sml.mat_. This file contains segmented nuclear pore complexes that you will be analyzing in the following steps.
:::{note}
You can check out all the pre-processing steps from _fitting raw data (raw camera frames)_ to _segmentation of NPCs_ in
**SMAP_manual_NPC.pdf** that can be downloaded [here](https://www.embl.de/download/ries/Documentation/).
:::

## Loading LocMoFit
1. Go to the **[ROIs]** tab.
2. Go to **[Evaluate]** tab and click **add module**.
3. In the popup window, select _LocMoFitGUI_ and click *ok*.

## Setup
We will be using an image model here ({doc}`more about model types<../basics/geometricModel>`).
1. Activate the _LocMoFitGUI_ module by clicking on it. Your SMAP window should look like this now:
   ![LocMoFit GUI in SMAP](../images/overview.png)
2. In the right panel, go to **[M1]** -> **[Model]**, click **load model**.
3. Navigate to the model directory, open _ring2d.png_. Now LocMoFit is ready to fit.

## Fitting
Fitting a site can be executed by clicking the site in the list of sites in the _ROI manager_ window.
1. Go to **[ROI]** -> **[Settings]**, click **show ROI manager**. This opens the **ROIManager** in a new window.
2. Go back to **[Evaluate]**, check **evaluate on** and **display** in the left panel.
	:::{Note}
	* The loaded evaluate plugins will be evaluated only when **evaluate on** is checked.
	* Result windows of the loaded evaluate plugins will be displayed only when **display** is checked.
	:::
3. In the _ROI Manager_ window, click on one site in the list of sites and wait for a few seconds. You should see a new window _LocMoFitGUI_, in which the localizations are plotted on the fitted model.
	:::{Note}
	The **LocMoFitGUI** window allows you to inspect the result right after fitting. The window may look different depending on the model type.
	:::
Now you have your first fit done! You can further explore a few sites to get familiar with the interface.

## Next tutorial
You are in the introductory series. The next tutorial is {doc}`Composite model<../tutorial/compositeModel>`.
