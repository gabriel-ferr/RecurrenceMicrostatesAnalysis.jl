#   RecurrenceMicrostatesAnalysis.jl for Devs
!!! tip
    All pull requests that introduce new functionality must be thoroughly tested and documented. Tests are required only for methods that you extend. We recommend reading the [Good Scientific Code Workshop](https://github.com/JuliaDynamics/GoodScientificCodeWorkshop).

    Always remember to add docstrings to your implementations, as well as tests to validate them.

##  RecurrenceMicrostatesAnalysis.jl backend
**RecurrenceMicrostatesAnalysis.jl** supports multiple backends, depending on the usage context. Each backend is implemented based on an [`RMACore`](@ref), which defines how the package computes a [`histogram`](@ref).

There are two main backend implementations:

- [`CPUCore`](@ref): defines how distributions are computed on the CPU. The default implementation is [`StandardCPUCore`](@ref).
- [`GPUCore`](@ref): defines how distributions are computed on the GPU. The default implementation is [`StandardGPUCore`](@ref).

```@docs
RMACore
CPUCore
GPUCore
StandardCPUCore
StandardGPUCore
```
!!! info
    Backend implementations are located in `src/core/cpu_core.jl` and `src/core/gpu/gpu_core.jl`. If you plan to implement a new backend, we recommend opening an [issue](https://github.com/DynamicsUFPR/RecurrenceMicrostatesAnalysis.jl/issues) on GitHub beforehand to discuss the design.

### Implementing an RMACore
Although it is possible to implement a custom [`RMACore`](@ref) directly, we **do not recommend** doing so. Instead, we strongly suggest implementing either a [`CPUCore`](@ref) or a [`GPUCore`](@ref).

This approach allows you to reuse utility functions such as `get_offsets` and `get_power_vector`, which expect an [`RMACore`](@ref) as input. Since these functions have different implementations for [`CPUCore`](@ref) and [`GPUCore`](@ref), writing a custom `RMACore` would require reimplementing them.

To avoid this, define a new struct that subtypes [`CPUCore`](@ref) or [`GPUCore`](@ref). In this case, the only required method to implement is [`histogram`](@ref). For example:
```@example mycore
using Random
using ComplexityMeasures
import RecurrenceMicrostatesAnalysis as rma
struct MyCore{M<:rma.MicrostateShape, S<:rma.SamplingMode} <: rma.CPUCore{M, S}
    shape::M
    sampling::S
end
```

```@example mycore
function histogram(
    core::MyCore,
    x::rma.StateSpaceSet,
    y::rma.StateSpaceSet
)
    
    # Construct the sampling space and determine the number of samples
    space = rma.SamplingSpace(core.shape, x, y)
    samples = rma.get_num_samples(core.sampling, space)

    # Precompute power vector and offsets
    pv = rma.get_power_vector(core, core.shape)
    offsets = rma.get_offsets(core, core.shape)

    # Allocate histogram
    hist = zeros(Int, rma.get_histogram_size(core.shape))

    # Task-local RNG (ignored for Full sampling)
    local_rng = TaskLocalRNG()

    # Histogram computation
    for m in 1:samples
        #   Get the sample.
        i, j = rma.get_sample(core, core.sampling, space, local_rng, m)
        #   Compute the microstate index.
        idx = rma.compute_motif(core.shape.expr, x, y, i, j, pv, offsets)
        @inbounds hist[idx] += 1
    end

    return Counts(hist, eachindex(hist))
end
```

!!! info
    To ensure compatibility with the internal API, custom backends must support the keyword argument `threads` (for [`CPUCore`](@ref)) or `groupsize` (for [`GPUCore`](@ref)), as required by the [`distribution`](@ref) overloads.

```@example mycore
data = rma.StateSpaceSet(rand(1000))
histogram(MyCore(rma.Rect(rma.Standard(0.27), 2), rma.SRandom(0.05)), data, data)
```

!!! warning
    CPU and GPU backends differ significantly in their execution models. In the GPU backend, random samples must be generated before histogram computation. The histogram itself is computed inside the `gpu_histogram!` kernel.

    See the [`StandardGPUCore`](@ref) implementation in `src/core/gpu_core.jl` for details.

!!! danger
    **RecurrenceMicrostatesAnalysis.jl** provides multiple backends that are only partially compatible.

    - CPU backends do not necessarily support spatial data.
    - Spatial analyses require dedicated implementations.
    - The GPU backend is fully incompatible with spatial data.

    Please consider these limitations carefully when extending or using backend functionality.

##  Adding a New Recurrence Function
### Steps
1. Define the mathematical expression of your recurrence function. It must return a binary value: `0` for non-recurrence and `1` for recurrence.
2. Define a new type `YourType <: `[`RecurrenceExpression`](@ref). Constant parameters (e.g., thresholds and metric) should be fields of this type.
3. Implement the appropriate [`recurrence`](@ref) dispatch:
   - Time series:  
     `recurrence(expr::YourType, x::StateSpaceSet, y::StateSpaceSet, i::Int, j::Int)`
   - Spatial data:  
     `recurrence(expr::YourType, x::AbstractArray{<:Real}, y::AbstractArray{<:Real}, i::NTuple{N,Int}, j::NTuple{M,Int})`
4. Add a docstring describing the mathematical definition and relevant references.
5. Add the recurrence expression to `docs/src/tutorial/recurrences.md`.
6. Add the expression to the [`RecurrenceExpression`](@ref) docstring.
7. Add tests to `test/distributions.jl` under the test set `recurrence expressions (with CPUCore)`.

!!! warning
    A recurrence function must always return `UInt(0)` or `UInt(1)`.

!!! todo
    To support GPU execution, also implement `gpu_recurrence(expr::YourType, x, y, i, j, n)`. See [`Standard`](@ref) for reference.

## Adding a New Sampling Mode

### Steps
1. Define how the sampling mode operates: which microstates are sampled, from which regions, and in what quantity. The [`SamplingSpace`](@ref) must be taken into account when designing the sampling logic.
2. Define a new struct that is a subtype of [`SamplingMode`](@ref). The struct may be empty (e.g. [`Full`](@ref)) or contain parameters such as a sampling rate (e.g. [`SRandom`](@ref)).
3. Implement the dispatch `get_num_samples(mode::YourType, space::SamplingSpace)` which determines the number of samples to be drawn given the sampling mode and the sampling space. Two sampling space types exist: `SSRect2` (time series) and `SSRectN` (spatial data).
4. Implement the dispatch `get_sample(core::RMACore, mode::YourType, space::SamplingSpace)` which returns the positions to be sampled. Separate implementations may be required for each [`RMACore`](@ref) and each [`SamplingSpace`](@ref). Full coverage is encouraged but not mandatory, provided that the supported cases are clearly documented in the docstring.
5. Add a docstring to your sampling mode describing its behavior and initialization. Follow the style of the existing sampling modes listed in [Implemented sampling modes](@ref).
6. Add your sampling mode to the list in `docs/src/tutorial/shapes_and_sampling.md`.
7. Add your type to the list in the [`SamplingMode`](@ref) docstring.
8. Add tests in `test/distributions.jl` under the test set `sampling mode (CPU backend)`.

!!! warning
    The `get_sample` logic differs between CPU and GPU backends. On the CPU, random samples are generated during histogram computation. On the GPU, samples must be generated beforehand, outside the kernel, and the kernel operates only on precomputed values.
    
##  Adding a new Microstate Shape
Defining a Microstate Shape is one of the most challenging tasks in this package (except for backend development, which is described by an [`RMACore`](@ref)).

A microstate shape acts as an intermediate structure between the sampling process and the recurrence function. Given an initial RP position $(i, j)$, it determines which additional recurrences must be evaluated and computes them using the recurrence expression. The resulting microstate is then converted into a decimal representation, which is used as an index in the histogram.

### Design considerations

Before implementing a microstate shape, it is essential to define its structure and reading order. For example, square microstates are typically read row-wise, while triangular microstates may be read column-wise. Each position in the microstate structure must be associated with a power of two in order to convert the binary microstate into a decimal index.
```math
\begin{pmatrix}
2^0 & 2^1 & 2^2 \\
2^3 & 2^4 & 2^5 \\
2^6 & 2^7 & 2^8
\end{pmatrix}
```

### Implementation steps
Define a new struct that is a subtype of [`MicrostateShape`](@ref). The struct must include a field `expr`, which stores the [`RecurrenceExpression`](@ref) used to compute recurrences at runtime.

Unlike [`RecurrenceExpression`](@ref) and [`SamplingMode`](@ref), a [`MicrostateShape`](@ref) does not require the implementation of `recurrence` or `get_sample` methods. Microstate computation is handled by the unified `compute_motif` function for the [`CPUCore`](@ref), and by `gpu_compute_motif` for the [`GPUCore`](@ref).

The only exception is spatial data, for which a custom `compute_motif` implementation is required. For example:
```julia
@inline function compute_motif(
    shape::RectN,
    x::AbstractArray{<: Real},
    y::AbstractArray{<: Real},
    idx::Vector{Int},
    itr::Vector{Int},
    power_vector::SVector{N, Int}
) where {N}
    
    index = 0
    dim = ndims(x) - 1
    copy!(itr, idx)

    @inbounds @fastmath for p in power_vector

        i = ntuple(k -> itr[k], dim)
        j = ntuple(k -> itr[dim + k], length(shape.structure) - dim)

        index += recurrence(shape.expr, x, y, i, j) * p

        itr[1] += 1
        for k in 1:length(shape.structure) - 1
            if (itr[k] > idx[k] + (shape.structure[k] - 1))
                itr[k] = idx[k]
                itr[k + 1] += 1
            else
                break
            end
        end
    end

    return index + 1
end
```

Although time-series microstate shapes do not require a custom `compute_motif` implementation, three utility functions must be defined to describe the properties of the shape:

1. `get_histogram_size(shape::YourType)`  
   Returns the length of the histogram, given by $2^\sigma$, where $\sigma$ is the number of recurrences in the microstate structure.

2. `get_power_vector(core::RMACore, shape::YourType)`  
   Returns the vector of powers of two used to convert the microstate into its decimal representation. This function differs between CPU and GPU backends due to integer size (`Int` on CPU, `Int32` on GPU).

3. `get_offsets(core::RMACore, shape::YourType)`  
   Returns the offsets relative to the initial position $(i, j)$ that define the remaining recurrence positions of the microstate. This function must be consistent with `get_power_vector` and also differs between CPU and GPU backends.

!!! tip
    For improved performance, we strongly recommend using `@generated` functions when implementing these utilities (except for spatial data).

Additionally, a [`SamplingSpace`](@ref) must be defined for the new microstate shape. For time-series data, this is typically `SSRect2`, while spatial data require `SSRectN`.

The sampling space must be initialized using the following constructor:
```julia
SamplingSpace(
    ::MicrostateShape, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{N, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{N, Float32}}}
)
```

After implementation, the following steps are required:

1. Add a docstring to your [`MicrostateShape`](@ref), explaining its behavior and initialization.
2. Add the definition to the section [Implemented microstates shapes](@ref) in `docs/src/tutorial/shapes_and_sampling.md`.
3. Add the type to the list in the [`MicrostateShape`](@ref) docstring.
4. Add tests to `test/distributions.jl` under the test set `motif shapes (with CPUCore)`.

### Example
First, we define the shape struct.
```@example line
using RecurrenceMicrostatesAnalysis
struct Line{N, B, E <: RecurrenceExpression} <: MicrostateShape
    expr::E
end

Line(expr::E, N; B = 2) where {E} = Line{N, B, E}(expr)
```

Next, we implements the three utils functions.
```@example line
@generated function get_histogram_size(::Line{N, B, E}) where {N, B, E}
    size = B^(N)
    return :( $size )
end
```
```@example line
@generated function get_power_vector(::CPUCore, ::Line{N, B, E}) where {N, B, E}
    expr = :(SVector{$N}( $([:(B^$i) for i in 0:(N-1)]... ) ))
    return expr
end
```
```@example line
@generated function get_offsets(::CPUCore, ::Line{N, B, E}) where {N, B, E}
    elems = [ :(SVector{2, Int}(0, $h)) for h in 0:(N - 1)]
    return :( SVector{$N, $(SVector{2, Int})}( $(elems...) ) )
end
```

Finally, we define our sampling space: (can be only for time series)
```@example line
using GPUArraysCore, StaticArrays
SamplingSpace(
    ::Line{N, B, E}, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{D, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{D, Float32}}}
) where {N, B, E<:RecurrenceExpression, D} = RecurrenceMicrostatesAnalysis.SSRect2(length(x), length(y) - N + 1)
```

And done, the shape can be used 😃. Remember to document it!

!!! warning
    The backend needs to have access to the util functions, them it is important to do the implementation inside the package module.

## Adding a New Quantifier
### Steps
1. Define a new quantifier type that is a subtype of [`QuantificationMeasure`](@ref).
2. Implement the corresponding [`measure`](@ref) dispatch used to compute the quantifier.
3. Add a docstring to the quantifier type, following the style of existing quantifiers.
4. Document the quantifier in `docs/src/tutorial/quantifiers.md`, including its definition, mathematical formulation, references, and examples when possible.
5. Add the quantifier to the list in the [`QuantificationMeasure`](@ref) docstring.
6. Add tests for the quantifier in `test/rqa.jl`.


## Adding a New GPU Metric
Since the [Distances.jl](https://github.com/JuliaStats/Distances.jl) package is not compatible with GPU execution, metric evaluations must be implemented manually to enable GPU support.

### Steps
1. Define a new type that is a subtype of [`GPUMetric`](@ref).
2. Implement the dispatch  
   `gpu_evaluate(::YourMetric, x, y, i, j, n)`  
   where `x` and `y` are `AbstractGPUVector{SVector{N, Float32}}`, `i` and `j` are indices, and `n` is the dimensionality of the vectors.
3. Add a docstring describing the metric, including its mathematical definition and parameters. Include references when applicable.
4. Document the metric in the section [Implemented GPU metrics](@ref) in `docs/src/tutorial/gpu.md`.
5. Add a reference to the metric in the [`GPUMetric`](@ref) docstring.

!!! danger
    GPU backends in **RecurrenceMicrostatesAnalysis.jl** operate exclusively with `Float32`. The use of `Float64` is not supported.
