#   Shapes and Sampling

This section describes the microstate shapes used in Recurrence Microstates Analysis (RMA)
and the sampling strategies employed to construct histograms and distributions.

##  Variations of microstates shapes
A microstate shape defines the local recurrence pattern extracted from a recurrence plot.
In **RecurrenceMicrostatesAnalysis.jl**, microstate shapes are represented by subtypes of
[`MicrostateShape`](@ref).

A `MicrostateShape` specifies:
- which relative positions are sampled to evaluate recurrences, and
- how the resulting binary recurrence pattern is mapped to a decimal representation.

The internal conversion from a microstate pattern to its decimal index is performed by the
`compute_motif` function.

```@docs
MicrostateShape
```

### Implemented microstates shapes
```@docs
Rect
Triangle
Diagonal
```

##  Sampling strategies
The sampling strategy determines which microstates are selected during histogram or
distribution construction.

Sampling behavior is defined by subtypes of [`SamplingMode`](@ref), while the set of valid
sampling positions is determined by a [`SamplingSpace`](@ref).

```@docs
SamplingMode
SamplingSpace
```

### Implemented sampling modes
```@docs
SRandom
Full
```