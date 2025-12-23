#   RecurrenceMicrostatesAnalysis.jl

```@docs
RecurrenceMicrostatesAnalysis
```

!!! todo "GitHub"
    RecurrenceMicrostatesAnalysis.jl is an open-source package available on [GitHub](https://github.com/DynamicsUFPR/RecurrenceMicrostatesAnalysis.jl).
    If you find this package useful, please consider giving it a star on GitHub and don't forget to cite [our work](https://doi.org/10.1063/5.0293708). 😉
    

##  About the documentation
RecurrenceMicrostatesAnalysis.jl documentation was designed to explain how to use the package while also introducing the
theoretical background of the RMA field. The bibliography used is listed in the [References](@ref)
section; **remember to cite the appropriate works if you use them**.

We begin this welcome section by introducing the [Input data for RecurrenceMicrostatesAnalysis.jl](@ref). It is very important
to understand the data types used by the package and their purposes before continuing with the rest of the documentation. We
also describe the [Output data from RecurrenceMicrostatesAnalysis.jl](@ref), explaining the type of data returned by the package when
computing distributions of recurrence microstates.

the **Tutorial** section explains how to use the package. We start with a brief introduction to the RMA framework and show how to
compute [Distributions](@ref) using RecurrenceMicrostatesAnalysis.jl. Next, we demonstate how to estimate RQA [Quantifiers](@ref)
using RMA and discuss several quantifiers defined specifically for RMA. This material forms the "basic level" of the documentation and is sufficient
to learn how to use the package effectively. We also include an introduction to [Operations](@ref) with recurrence microstates.

If you want to learn more about RecurrenceMicrostatesAnalysis.jl, the [Recurrence Functions](@ref) section discusses variations in
computing recurrence between two states, and the [Shapes and Sampling](@ref) section provides explanations about different motif shapes, which
are used to extract specific information from an abstract Recurrence Plot.

We also provide a pipeline for [GPU](@ref) computations, which we recommend reading if you intend to use this framework. Moreover, the
[Performance Tips](@ref) section offers advice on improving the performance of RecurrenceMicrostatesAnalysis.jl and avoiding
common pitfalls.

This documentation also includes some sections with applied examples:
- [RMA with Machine Learning](@ref)

Finally, if you are a developer interested in contributing to RecurrenceMicrostatesAnalysis.jl, we recommend reading the section
[RecurrenceMicrostatesAnalysis.jl for Devs](@ref).

##  Input data for RecurrenceMicrostatesAnalysis.jl

RecurrenceMicrostatesAnalysis.jl accepts two types of input (each of them with a different backend):

- [`StateSpaceSet`](@ref) — used for multivariate time series, datasets, or state-space sets. It is employed by the backend when working with Recurrence Plots (RP) or Cross-Recurrence Plots (CRP). If you are working with RP or CRP, we strongly recommend using this data type, since the backend is optimized for it in this context.

- `AbstractArray{<: Real}` — used for spatial data. This allows RMA to be applied in the generalized framework of Spatial Recurrence Plots (SRP) [Marwan2007Spatial](@cite). If you provide a `Matrix`, this input type can also be used instead of a [`StateSpaceSet`](@ref); however, we do **not** recommend it, as the backend for `AbstractArray{<: Real}` is heavier and incompatible with some features.

- `AbstractGPUVector` - used for analysis of time series using the GPU backend. We explain how it works better in the section [GPU](@ref).

!!! warning
    RMA with SRP is an open research field. We include this functionality in the package for exploratory purposes, but the method is not
    yet mature enough for production use. Nevertheless, feel free to experiment with it in your research. 😃

```@docs
StateSpaceSet
```

##  Output data from RecurrenceMicrostatesAnalysis.jl

When computing the RMA distribution, RecurrenceMicrostatesAnalysis.jl returns a [`Probabilities`](@ref) structure. This type is
provided by [ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/), allowing this package
to interoperate naturally with its tools and workflows.

```@docs
Probabilities
Counts
```