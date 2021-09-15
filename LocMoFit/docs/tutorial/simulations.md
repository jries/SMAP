# Simulating SMLM data

In this tutorial you will learn how to simulate SMLM data based on a geometric model. You will need:
* **SMAP** installed. Further information can be found on our [GitHub](https://github.com/jries/SMAP/) site.

## Preparation
You have to load the geometric model for generating localizations. Here we use the point model of the nuclear pore complex as an example.
1. Start SMAP.
2. Load a new instance of the plugin *LocMoFitGUI* (see {doc}`quick start<../quickstart>` to recap).
{doc}`quick start<./basics/quickstart>`
{doc}`quick start<../basics/quickstart>`
3. Then we load the individual models and set up the arguments of the model parameters:
	* Go to **[M1]** -> **[Model]**, click **load model**. Navigate to the default model directory *'SMAP/LocMoFit/models'* (*SMAP* is the root folder of SMAP), open _NPCPointModel_flexible2.m_.
	
## Setup
Once the model is loaded, it is ready for the simulation engine.
1. We first connect the engine to the model:
	* Go to [ROIs]->[Segment], click *SimulateSites* in the list of plugin.
	* Click **Connect to LocMoFit** and check LocMoFitGUI in the new window.
2. Next, you have to define your desired values of model parameters.
	* Click **Set model pars**.
	* In the new window, define the value for a parameter (can be identified in the columns *Name*, *Type*, and *Model*) in the column *Value*.
	* Newly inputted values is automatically saved.
	:::{note}
	Instead of a constant value, a parameter can be a random variable of a uniform distribution with specified boundaries. To enable this, you just have to input two values, separated by a space, in the column *value*. For example, *'-15 15'* defines the lower and upper boundaries as -15 and 15 respectively.
	:::
	:::{note}
	You can also associate model parameters. See {doc}`this note<../applicationNote/simulation>` for more information.
	:::
3. Then you have to set up the SMLM properties.
	* Go back to the GUI of *SimulateSites*.
	* Input the desired values for the SMLM properties. You can hover over the respective properties to get their definitions.
## Simulation
Click **Run** in the GUI of *SimulateSites*. The simulation is done when *'ROIManager.Segment.SimulateSites finished'* shows up in the status bar. This usually takes a few seconds.
