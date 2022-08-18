Geometric model classes
=======================
These most basic classes allow defining a geometric model in user-defined forms (see :class:`geometricModel<@geometricModel.geometricModel>`) or in a parametric form (see :class:`parametricModel<@parametricModel.parametricModel>`). 

Geometric model
---------------
.. important::

	All the direct subclasses of :class:`geometricModel<@geometricModel.geometricModel>` are required to contain the method :meth:`reference`.

.. method:: reference(obj, par, dx)
	
	This funtion samples coordinates of the model as reference.
	
	Usage:
		[model, p]= reference(obj, par, dx)
		
	Input:
	  * **obj** (:class:`geometricModel<@geometricModel.geometricModel>`) – an object of any subclass of :class:`geometricModel<@geometricModel.geometricModel>`.
	  * **par** (structure array) – each field contains a parameter value and its fieldname should be the parameter name.
	  * **dx** (numeric scalar) – sampling rate.

	Output:
	  * **model** (structure array) – a structure object. Its fieldnames are x, y, z, and n, indicating the amplitudes n at xyz positions of the sampled model points.
	  * **p** (structure array) – additional information of the model.

.. method:: getDerivedPars(obj, pars)

	This funtion calculates additional parameters derived from the geometric paramaters.

	Usage:
		derivedPars = getDerivedPars(obj, pars)
		
	Input:
	  * **obj** (:class:`geometricModel<@geometricModel.geometricModel>`) – an object of any subclass of :class:`geometricModel<@geometricModel.geometricModel>`.
	  * **pars** (structure array) – each field contains a parameter value and its fieldname should be the name of the geometric parameter.

	Output:
	  * **derivedPars** (structure array) – each field contains a parameter value and its fieldname should be the name of the derived parameter.
	 
.. automodule:: @geometricModel

.. autoclass:: geometricModel
    :members:

Parametric Model
----------------
.. automodule:: @parametricModel

.. autoclass:: parametricModel
    :members: