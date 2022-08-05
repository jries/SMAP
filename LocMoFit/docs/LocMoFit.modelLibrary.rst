Model library
=============
All the names of the models follow the rule: **'[geometry][dimension][p]_[parameterization]_[...]_[model form]'**, with the simplest case as '[geometry][dimension]' (e.g., :class:`ring3D<models.ring3D>`). For example, :class:`dualEllipse3D_avgR_discrete<models.dualEllipse3D_avgR_discrete>` means a dual-ellipse geometry in 3D, parametrized by the average radius, in a discrete form.

.. note::
	The individual components of the names

	.. hlist::
		:columns: 1
		
		* [geometry]: the geometry.
		* [dimension]: the model dimension.
		* [p]: (optional) when it is mentioned, the model is parametric, otherwise not. 
		* [parameterization]: (optional) when multiple parameterizations exist, this is added to identify the specific parameterization. When the [parameterization] is not mentioned, the model is assumed to be with the default parameterization. 
		* [...]: (optional) if the [parameterization] itself is not sufficient to differentiate the implementations, extra labels are added here.
		* [model form]: (optional) the form of how the model is implemented. When it is not provided, the model is in the continuous form.

.. automodule:: models

Arc
---
.. autoclass:: arc2D
	:members:
	:show-inheritance:

.. autoclass:: arc2D_arcLen
	:members:
	:show-inheritance:

Bucket
------
.. autoclass:: bucket2D
	:members:
	:show-inheritance:

Ring
----
.. autoclass:: ring2D
	:members:
	:show-inheritance:

.. autoclass:: ring3D
	:members:
	:show-inheritance:

Spline
------
.. autoclass:: cspline3D_midPoint
	:members:
	:show-inheritance:

Tube and derivitives
--------------------
.. autoclass:: csplineClosedTube3D_midPoint
	:members:
	:show-inheritance:

.. autoclass:: csplineTube3D_midPoint
	:members:
	:show-inheritance:

.. autoclass:: csplineTube3D_xyz
	:members:
	:show-inheritance:

Dual rings and derivitives
--------------------------
.. autoclass:: dualEllipse3D_avgR_discrete
	:members:
	:show-inheritance:

.. autoclass:: dualEllipse3D_discrete
	:members:
	:show-inheritance:

.. autoclass:: dualRing3D_discrete
	:members:
	:show-inheritance:

Sperical model and derivitives
------------------------------
.. autoclass:: sphericalCap3D_surfaceArea
	:members:
	:show-inheritance:

.. autoclass:: sphericalCap3Dp_surfaceArea
	:members:
	:show-inheritance:

.. autoclass:: spheroid3Dp_surfaceArea
	:members:
	:show-inheritance:

.. autoclass:: spheroidCap3Dp_surfaceArea
	:members:
	:show-inheritance:
	
2D projection of 3D geometry
----------------------------
.. autoclass:: hemispheroid2D
	:members:
	:show-inheritance:

.. autoclass:: thickRing2D
	:members:
	:show-inheritance:

Random geometry
---------------
.. autoclass:: locsBG3D
	:members:
	:show-inheritance:
