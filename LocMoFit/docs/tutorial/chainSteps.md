# Chaining steps

In this tutorial you will learn how to chain fitting steps with the GUI. We will set up three fitting steps with different models, using the results of the previous step as the starting parameter for the next one. You will need:
* **SMAP** installed
* _U2OS_Nup96_BG-AF647_demo_sml.mat_
* _NPC3D_step2_ring_SMLMModelFit.mat_
* _NPC3D_step3_points_SMLMModelFit.mat_

:::{important}
Please first finish the tutorial {doc}`composite model<./compositeModel>`.
:::

## Loading LocMoFit
In the previous tutorial you have loaded one **SMLMModelFitGUI** plugin. We need two more now:
* Go to **[Evaluate]** tab and click on **add module**.
* In the popup window, select _SMLMModelFitGUI_ and click **ok**.
* Repeat the two steps above once more.

## Setup
Next, we load the two following steps:
1. In the left panel, click *SMLMModelFitGUI_2*.
2. On the right panel, go to **[Settings]**, click **load**, navigate to the settings directory, and select *NPC3D_step2_ring_SMLMModelFit.mat*.
3. In the left panel, click *SMLMModelFitGUI_3*.
4. On the right panel, go to **[Settings]**, click **load**, navigate to the settings directory, and select *NPC3D_step3_points_SMLMModelFit.mat*.
:::{note}
The different steps are connected via the parameter converter in the **[Convert]** tab. It calculates values according to the **rules** based on its **source** and then writes the values to the **target_fit** (see the column names of the convert table).
:::

## Fitting
Click on one site in the _ROI manager_ window. Now you should see three view windows, one for each step.

## Batch fitting
After inspecting several sites, we are now moving on to fit all sites. Such a batch analysis can be executed with the _redraw all_ function:
Go to the **[Evaluate]** tab. In the left panel, uncheck **display** and click **redraw all**. You will see the analysis going down the list of sites.
:::{Note}
* While running, you should see the message _"redrawall: site [current site] of 100"_ in the status bar.
* You will know the analysis is done when _"redrawall: completed"_ shows up in the status bar.
:::

## Summary of parameter values
To get a summary of parameter values, we can use the SMAP plugin _summarizeModFitNPC3D_:
* Go to the drop-down menu **[Plugins]** -> **[ROIManager]** -> **[Analyze]** -> **[_summarizeModFitNPC3D_]**
* Click **Run** in the new window.

## Next tutorial
This is the end of the introductory series. If you started from {doc}`quick start<./quickstart>` and followed along the _Next tutorial_ section in each tutorial, you should be ready for more advanced analyses with LocMoFit.
