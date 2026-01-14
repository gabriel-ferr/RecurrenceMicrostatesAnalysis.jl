#   CHANGELOG

This changelog is maintained starting from version 0.4.0.

##  0.4.0
- Complete redesign of the backend using abstract types, structs, and function overloading, replacing the previous function-based implementation.
- The input type of the `distribution` function has been changed to `StateSpaceSet` for time series, improving interoperability with [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/).
- The output type of the `distributions` function has been changed to `Probabilities`, improving interoperability with [ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/)
- Added a new GPU backend based on [KernelAbstractions.jl](https://juliagpu.github.io/KernelAbstractions.jl/stable/).
- Implemented operations for computing microstate permutations.

### BREAKING CHANGES
Version 0.4.0 is not backward compatible with v0.3.1. The backend has been redesigned from a function-based architecture to one based on abstract types, structs, and method overloading.

The overall usage of the `distribution` function remains similar, but its configuration options have changed. Significant changes were also introduced in RQA computation: the package now uses a unified `measure` interface and the `QuantificationMeasure` type instead of the standalone functions used previously.

**Key changes:**
- The functions `rentropy`, `rrate`, `determinism`, `laminarity`, and `disorder` have been replaced by the generic `measure` function together with the corresponding subtypes of `QuantificationMeasure`.
- The function `find_parameters` has been replaced by `optimize` using the `Threshold` type.
- The function `prepare` has been removed, since this functionality is now provided by [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/).
- The `distribution` function is no longer configured via **keyword arguments** for the microstate shape and recurrence function. These are now specified using the types `MicrostateShape` and `RecurrenceExpression`. 
    The sampling mode and sampling ratio remain configurable via **keyword arguments**, but `SamplingMode` is now a type rather than a symbol, and the keyword `sampling_mode` has been renamed to `sampling`.
- The inputs of `distributions` and `histogram` for time series are no longer `Vector` or `Matrix`. They must now be provided as a `StateSpaceSet`.