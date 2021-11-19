# Parameter table
An example parameter table:

![m1_Parameters](../images/m1_Parameters.PNG)

## Column
* _name_: parameter names.
* _value_: initial values of parameters.
* _fix_: specify the parameter should be fix or or not. If checked, the parameter will be set to the value you defined and will be a constant that are not fitted.
* _lb_: relative lower boundaries of parameter ranges.
* _ub_: relative upper boundaries of parameter ranges.
* _type_: types of parameters. This is not editable.
* _min_: absolute lower boundaries of parameter ranges.
* _max_: absolute upper boundaries of parameter ranges.
* _label_: user-defined labels for the corresponding parameters.

   :::{note}
   For example, _value_ 30, _lb_ -10, _ub_ 20, _min_ 10, and _max_ 45 results to a parameter range of [20 45]. This is based on first get [value+lb value+ub] = [20 50], and then check whether this is beyond [min max] = [10 45] or not. If this is the case, the range will be set to the min or max values so the final upper boundary is 45 but not 50.
   :::