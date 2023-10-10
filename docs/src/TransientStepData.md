# TransientStepData


## *What is it?*

TransientStepData objects hold information pertaining to *transient sorption steps*. This is dimensionless data that describes the transition from one equilibrium to another (e.g., concentration, dilation, etc) as a function of time. 


## *How do I make one?*

```@docs
TransientStepData
```

A number of functions exist to make your life easier with handling this data. 

```@docs
resample(::TransientStepData, ::Any, ::Any)
dataset(transient::TransientStepData)
strip_measurement_to_value(meas::TransientStepData)
```