# Multicolor protein complex reconstruction

## Overview

In this tutorial, you will learn how to perform the reconstruction of a multi-color protein complex based on multiple dual-color datasets.

## Procedure

### Loading data
On your PC, you can find the project folder that TAs opened for you.
1. **Copy** the path to this folder.
2. Open MATLAB and locate to the SMAP folder (_you will be in this folder by default_).
3. In the MATLAB console, type **SMAP** and hit the **enter** key.
4. Wait until you see _'all initialized'_ in the **status bar**.
5. Click on **Load**, go to the folder path you copied. Go to the subfolder **'data'**, select the file named after your group.

### Assign colors to localizations
1. Assign the localizations to the different color channels:
	* Got to **[Process]** -> **[Assign2C]** -> **Intensity2ManyChannel**.
	* In the fields, select _phot1_ and _phot1_ respectively.
	* Click on ROI1. In the point cloud, box out the population on the left.
	* Click on ROI2. In the point cloud, box out the population on the right.
	* Adjust the boxes to ensure good separation.
	* Click **Run**.
2. To display the 2 different channels, go to the **[Render]** tab and add a layer by clicking on **[+]** Change the Channel (**Ch**) to 2 and change the LUT to a different one than used in Layer 1
and click on **Reconstruct**.
3. Filter the localizations in **[Layer1]** according to TAs' instructions.
4. Uncheck the _layer 2_.

### Segmentation of NPCs
The filtered and drift corrected superresolution image can be segmented for NPC automatically.
1. Automatic segmentation
	1. Open the _ROI manager_.
	2. **Run** the _makeCellGrid_ plugin in **[Segment]**.
	3. **Redraw all** in **[Settings]**.
	4. Run the _segmentNPC_ plugin in **[Segment]**. Typical parameters for NPC
	segmentation are cutoff: 0.06, diameterNPC: 100, rim: 20. Check the box saveon
	before running the plugin. If only a part of the image should be segmented,
	additionally check getmask.
	5. In **[Evaluate]** check only _NPCsegmentCleanup_. For _NPCsegmentCleanup_ standard
	parameters are suitable in most cases. However, max average PSF can be set to
	140 to further exclude out of focus localizations.
	6. Run the plugin by checking _evaluate on_ (_display_ can be unchecked to increase
	processing speed) and a click on **redraw all**. This step can take several minutes.
	7. Switch to the ROI manager and click **update** above the list of ROIs. ROIs which were
	unsuccessful during the _NPCsegmentCleanup_ quality check will be shown with a – at
	their end.
	8. In **[Helper]** set the 1st sort to _descend_. Following that, choose Annotation in the dropdown menu. A list opens upon clicking **select**. Set **list4->value** and start with **Run**.
	This step sorts the list of ROIs in the ROI manager based upon the results of the
	quality check through the _NPCsegmentCleanup_ plugin.
	9. In the ROI manger, select all ROI names with a – at their end and delete them by
	clicking **Remove** or pressing the delete key on your keyboard.

### Estimation of spatial parameter
Next we fit a geometric model to single sites. This process searches in the parameter space of the model to match the model to the data. The parameter values that yields the best match are the estimates of the parameters.
1. Loading the plugin **SMLMModelFitGUI**:
	* Go to **[Evaluate]** tab and click on *add module*.
	* In the popup window, select _SMLMModelFitGUI_ and click *ok*.
	* Repeat the previous 2 steps 2 times.
2. Loading settings for model fitting:
	* Load corresponding settings to each plugin in the following table:

	| Plugin | Settings |
	| :------- | :---- |
	| SMLMModelFitGUI | Step1_SMLMModelFit.mat |
	| SMLMModelFitGUI_2 | Step2_SMLMModelFit.mat |
	| SMLMModelFitGUI_3 | Step3_SMLMModelFit.mat |

	* In the left panel, **click** on the corresponding plugin.
	* In the right panel, go to the **[Settings]** tab.
	* Click on **load**. Select the corresponding settings in the sub-folder _settings_. Click on **open**.
	* Repeat the previous 3 steps for the rest of plugins in the table.

3. Perform fitting on one site at a time:
	* Go to _ROI manager_ window.
	* Click any site in the list of sites.
	* Inspect the result in the popup windows.
4. Perform fitting through all sites:
	* Go to **[Evaluate]** tab. In the left panel, uncheck **display** and click **redraw all**.
	* Wait until all the first 200 sites are analyzed.
	* Click on **Stop** in the main SMAP window.
5. In MATLAB, open the script **summary** from the _script_ folder, run the script, and note down the _azimuthal shift_ displayed in the console.
6. In **[Helper]** set the 1st sort to _ascend_. Following that, choose Annotation in the dropdown menu. A list opens upon clicking **select**. Set **list3->value** and start with **Run**.
7. In the _ROI Manager_, delete all sites with the value 8 or 9 in the third digit of the number at the end of their site names.

### Generation of the average NPC per dataset
1. Use another SMLMModelFitGUI entity to perform averaging:
	* Go to **[Evaluate]** tab and click on **add module**.
	* In the popup window, select **SMLMModelFitGUI** and click on ok.
	* Click on *SMLMModelFitGUI_4* in the left panel.
	* Go to the **[Settings]** tab in the right panel and click **load**.
	* In the popup window, open the **Step3_SMLMModelFit.mat** in the folder **/settings**.
	* Go to **[convert]** tab in the right panel and click on **Match**.
	* In the first row, change **this step** to _SMLMModelFitGUI_3_.
	* In the list bellow, uncheck _m1.zrot_ and _m1.azimuthalShift_ in the popup window. Click on **Apply**.
	* Go to **[M1]** tab in the right panel. Fill the _azimuthal shift_ you noted down in the _value_ column of _azimuthalShift_.
	* Go to **settings** tab in the right panel, check **Run alignment**
	* Click on **redraw all** in the left panel.
	* Wait until all the sites are analyzed.
2. Register all sites to the average:
	* In the dropdown menu, open **[Process]** -> **[Modify]** -> **DefineMainCoordinates**
	* Check **xnm**, **ynm**, **znm**, and **field**.
	* In the first row, select **xnmaligned_SMLMModelFitGUI_4**, **ynmaligned_SMLMModelFitGUI_4**, and **znmaligned_SMLMModelFitGUI_4** from left to right.
	* Click on **RUN**, and wait until **Process.Modify.DefineMainCoordinates finished** show up in the status bar.
	* Double right click on the window **[Reconstructed superresolution image]**, you will see the average.
3. Save the average:
	* Remove all filtering: uncheck **[filt]** for each row in the **localization table**.
	* Go to **[File]** tab and check **only save visible**.
	* Click on **Save** and type the file name in the popup window and **Save**.

### Multi-color reconstruction.
1. Open a new SMAP session.
2. Load all dual-color averages.
3. Define different layers
