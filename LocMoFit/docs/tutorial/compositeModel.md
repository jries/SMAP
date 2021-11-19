# Composite model

In this tutorial you will learn how to build a composite model with the GUI. We will build a 3D dual-ring model by combining two times the identical ring model. You will need:
* **SMAP** installed. Further information can be found on our [GitHub](https://github.com/jries/SMAP/) site.
* _U2OS_Nup96_BG-AF647_demo_sml.mat_
* _ring3d.mat_
* _dualRing_model1_fitPar.csv_
* _dualRing_model2_fitPar.csv_

The four files above can be downloaded [here](https://www.embl.de/download/ries/LocMoFit/).

## Preparation
:::{important}
Skip step 1 and 2 if you continue from {doc}`quick start<../basics/quickstart>`.

Instead, remove the previous loaded **LocMoFitGUI** if you continue from {doc}`quick start<../basics/quickstart>`:
	1. Go to **[ROIs]** -> **[Evaluate]** and click on the module to remove.
	2. Click **remove**.
:::
1. Start SMAP.
2. Load the dataset _U2OS_Nup96_BG-AF647_demo_sml.mat_.
3. Load a new instance of the plugin **LocMoFitGUI** (see {doc}`quick start<../basics/quickstart>` if you forget how to do it).

## Setup
We will combine two identical rings (in the image format) in 3D to form a 3D dual-ring model ({doc}`more about model types<../basics/geometricModel>`).

1. We first load the individual models and set up the arguments of the model parameters:
	* Go to **[M1]** -> **[Model]**, click **load model**. Navigate to the model directory, open _ring3d_img.mat_. Go to the tab **[Parameters]** and click the button **load**. In the new window, navigate to the settings directory and open dualRing_model1_fitPar.csv. Click **Import** button in the popup window.
	* Click **[+]** next to **[M1]** to add a second model tab **[M2]**. In **[M2]**, click **load model**. Navigate to the model directory, open _ring3d_img.mat_ again. Go to the tab **[Parameters]** and click the button **load**. In the new window, navigate to the settings directory and open dualRing_model2_fitPar.csv. Click **Import** button in the popup window.
	:::{note}
	In this tutorial we loaded the pre-defined parameter arguments to simplify the procedure. For all the structures we fit in our manuscript, we provide optimized arguments. However, when you work on a new structure, you have to tweak the arguments yourself and modify arguments in the **parameter table** (find out more {doc}`here<../basics/parTable>`).
	:::
2. We then set up the **optimizer**:
	* In the right panel, go to **[Settings]** (next to **[M1]**), select _particleswarm_ for the **optimizer**.
	* Click the **+** button below the (empty) parameter list to add a new row in the parameters table of the optimizer. Select _UseVectorized_ in the column **parameters** and fill in _true_ in the column **Value**.
3. To estimate the starting parameters, we use **Convert**, which converts a rule to the starting value of a certain parameter.
	Go to tab **[Convert]** and fill in the table as the following (you can use the **+** button to add a new row):
	| Source    | Rule             | Target\_fit | Target\_usr |
	| --------- | ---------------- | ----------- | ----------- |
	| this step | median(locs.xnm) | m1.lPar.x   |             |
	| this step | median(locs.ynm) | m1.lPar.y   |             |
	| this step | median(locs.znm)+40 | m1.lPar.z   |             |
	:::{note}
	Here we assign the median (xnm, ynm, and znm refer to the values on the respective coordinate axis, in nanometer) position of localizations (locs) as the initial guess for the center position of the model (e.g., m1.lPar.x means the x position of model 1).
	:::

4. To inspect the model with the starting parameters, we use the functions **pick site** and **preview**.
	* Go to **[M1]** and click **pick site**.
	:::{note}
	**pick site** allows you to click on a site without executing fitting.
	:::
	* Click on one site in the list of sites in the _ROI manager_ window.
	* Click **preview**. A view window should pop up. You can see different views in tabs **[xy]**, **[xz]** and **[yz]**.

## Fitting
Click on the same site in the _ROI manager_ window. You should see the view window updated.

## Next tutorial
You are in the introductory series. The next tutorial is {doc}`Chaining steps<./chainSteps>`
