# Geometric model
:::{note}
This page is still incomplete. We are working on the content now and will update it soon.
:::
## Model format
* Image: _png_ and _mat_ files.
* Matlab function: _m_ files.

## Model type
* Image: a rigid template with only extrinsic parameters.
* Continuous: a model with both intrinsic and extrinsic parameters. During the optimization, a continuous model generates a density map and the likelihood value of each localization is calculated based on spatial interpolation on the map. Since the map is the same for all localizations, they share the same uncertainty (usually its mean values).
* Discrete: a model with both intrinsic and extrinsic parameters. During the optimization, a discrete model generates a set of model points and the likelihood value of each localization is calculated based on its distance to each model point. The real uncertainty of each localization is used.
