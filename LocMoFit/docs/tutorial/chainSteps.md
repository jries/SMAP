# Chaining steps

:::{note}
Time required: ~10 min.
:::

:::{note}
In this tutorial we loaded the pre-defined parameter arguments to simplify the procedure. For all the structures we fit in our manuscript, we provide optimized arguments. However, when you work on a new structure, you have to tweak the arguments yourself and modify arguments in the **parameter table** (find out more {doc}`here<../basics/parTable>`).
:::
	
## Task
Chaining fitting steps with the GUI. We will set up two fitting steps with different models, using the results of the previous step as the initial parameters for the next one.

## Requirment
* Software: **SMAP** installed. Further information can be found on our [GitHub](https://github.com/jries/SMAP/) site.
* Localization data: _U2OS_Nup96_BG-AF647_demo_sml.mat_
* Fitting settings:
	* _NPC3D_step1_ring_LocMoFit.mat_ (the one you saved in the last [tutorial](./compositeModel#Saving_all_settings))
	* _NPC3D_step2_points_LocMoFit.mat_

The other two files can be downloaded [here](https://www.embl.de/download/ries/LocMoFit/).

:::{important}
Please first finish the tutorial {doc}`composite model<./compositeModel>`.
:::

## Preparation
1. Start SMAP ({doc}`how to? <../howto/SMAP.runSMAP>`).
	:::{important}
	If you continue from the previous tutorial, please close the old one and start a new session.
	:::
2. Load the dataset _U2OS_Nup96_BG-AF647_demo_sml.mat_.
3. Load an instance of the plugin **LocMoFitGUI** (see {doc}`quick start<../tutorial/quickstart>` if you forget how to do it).

## Loading LocMoFit
Now you need two instances of LocMoFit:
* Go to **[Evaluate]** tab and click on **add module**.
* In the popup window, select _LocMoFitGUI_ and click **ok**.
* Repeat the two steps above once more.

Now you should see the **LocMoFitGUI** and **LocMoFitGUI_2** in the loaded moduals. These are for two different fitting steps respectively.

## Setup
Next, we load the LocMoFit settings for fitting two different models (step 1: {class}`ring3D<models.ring3D>`; step 2: {class}`dualRing3D_discrete<models.dualRing3D_discrete>`):
1. In the left panel, click *LocMoFitGUI*.
2. On the right panel, go to **[Settings]**, click **load**, navigate to the settings directory, and select _NPC3D_step1_ring_LocMoFit.mat_ (which you saved earlier).
3. In the left panel, click *LocMoFitGUI_2*.
4. On the right panel, go to **[Settings]**, click **load**, navigate to the settings directory, and select _NPC3D_step2_points_LocMoFit.mat_.
5. Go to tab **[Convert]**, you will find the table is pre-defined to simplify the tutoral. To better explain what's in the table, we left one row for you to fill in:
	* click **+** button to add a new row, and then fill in:
	
	| Source    | Rule             | Target\_fit | Target\_usr |
	| --------- | ---------------- | ----------- | ----------- |
	| LocMoFitGUI | pars.m2.lPar.z | m1.mPar.ringDistance | 	 |
:::{note}
**[Convert]** can be used to convert the fitted values in the previous step to a initial parameter of the current step. It calculates values according to the **rules** based on its **source** and then writes the values to the **target_fit** (see the column names of the convert table). For example, you just defined to assign the parameter m1.mPar.ringDistance in this step from the **Source** _LocMoFitGUI_ based on the value calculated by the **Rule** pars.m2.lPar.z.

_m1.mPar.ringDistance_ and _pars.m2.lPar.z_ are IDs of the corresponding parameters. m1.mPar.ringDistance means the parameter _ringDistance_, which is a intrinsic parameter (_mPar_) of model 1 (_m1_). See the syntax that can be used in **[convert]**.
:::

## Fitting
Click on one site in the _ROI manager_ window. Now you should see two view windows, one for each step.

## Batch fitting
After inspecting several sites, we are now moving on to fit all sites. Such a batch analysis can be executed with the _redraw all_ function:
Go to the **[Evaluate]** tab. In the left panel, uncheck **display** and click **redraw all**. You will see the analysis going down the list of sites.
:::{Note}
* While running, you should see the message _"redrawall: site [current site] of 100"_ in the status bar.
* You will know the analysis is done when _"redrawall: completed"_ shows up in the status bar. This usually takes around 1-2 minutes.
:::

## Summary of parameter values
To get a summary of parameter values, we can use the SMAP plugin _summarizeModFitNPC3D_:
* Go to the drop-down menu **[Plugins]** -> **[ROIManager]** -> **[Analyze]** -> **[_summarizeModFitNPC3D_]**
* Click **Run** in the new window.

## The end of the introductory series
This tutorial is the end of the introductory series. If you started from {doc}`quick start<./quickstart>` and followed along the _Next tutorial_ section in each tutorial, you should be ready for more advanced analyses with LocMoFit.

## Further reads
* 
*