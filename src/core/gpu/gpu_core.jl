export GPUCore, StandardGPUCore

##########################################################################################
#   RMACore: GPU
##########################################################################################
"""
    GPUCore{B, M<:MicrostateShape, S<:SamplingMode} <: RMACore

Abstract GPU backend that implements the **RecurrenceMicrostatesAnalysis.jl** execution pipeline on
graphics processing units.

The package provides a default implementation via [`StandardGPUCore`](@ref).

Concrete subtypes of `GPUCore` must define the following fields:
- `backend`: the GPU backend device (e.g. `CUDABackend`, `MetalBackend`).
- `shape`: the [`MicrostateShape`](@ref) used to construct microstates.
- `sampling`: the [`SamplingMode`](@ref) used to sample the recurrence space.

#   Implementations
- [`StandardGPUCore`](@ref)
"""
abstract type GPUCore{B, M<:MicrostateShape, S<:SamplingMode} <: RMACore end

##########################################################################################
"""
    StandardGPUCore{B, M<:MicrostateShape, S<:SamplingMode} <: GPUCore{B, M, S}

Default GPU backend implementation for **RecurrenceMicrostatesAnalysis.jl**.

This type provides the standard execution pipeline for computing recurrence microstate distributions
on GPU devices.

#   Initialization
```julia
core = GPUCore(backend, shape, sampling)
```
"""
struct StandardGPUCore{B, M<:MicrostateShape, S<:SamplingMode} <: GPUCore{B, M, S}
    backend::B
    shape::M
    sampling::S
end

##########################################################################################
GPUCore(backend::B, shape::M, sampling::S) where {B, M<:MicrostateShape, S<:SamplingMode} = StandardGPUCore(backend, shape, sampling)

##########################################################################################
#   Implementation: compute_motif
##########################################################################################
@inline function gpu_compute_motif(expr, x, y, i, j, power_vector, offset, n)
    index = zero(Int32)

    @inbounds begin
        for m in eachindex(power_vector)
            dw, dh = offset[m]
            @fastmath index += power_vector[m] * gpu_recurrence(expr, x, y, i + dw, j + dh, n)
        end
    end

    return @fastmath index + 1
end

##########################################################################################
#   Implementation: histogram
##########################################################################################
#   Based on time series: (GPU)
#.........................................................................................
"""
    histogram(core::StandardGPUCore, [x], [y]; kwargs...)

Compute the histogram of recurrence microstates for an abstract recurrence structure constructed
from the input data `[x]` and `[y]`.

If `[x]` and `[y]` are identical, the result corresponds to a Recurrence Plot (RP); otherwise, it
corresponds to a Cross-Recurrence Plot (CRP).

The result is returned as a [`Counts`](@ref) object representing the histogram of recurrence
microstates for the given input data.

!!! note
    The resulting histogram is copied from GPU memory back to the CPU.

This method implements the GPU backend using a [`GPUCore`](@ref), specifically the
[`StandardGPUCore`](@ref) implementation.

### Arguments
- `core`: A [`StandardGPUCore`](@ref) defining the GPU backend configuration.
- `[x]`: Input data provided as an `AbstractGPUVector`.
- `[y]`: Input data provided as an `AbstractGPUVector`.

### Keyword Arguments
- `groupsize`: Number of threads per GPU workgroup.

### Examples
```julia
using CUDA
gpudata = StateSpaceSet(Float32.(data)) |> CuVector
core = GPUCore(CUDABackend(), Rect(Standard(0.27f0; metric = GPUEuclidean()), 2), SRandom(0.05))
dist = histogram(core, gpudata, gpudata)
```
"""
function histogram(
    core::StandardGPUCore,
    x::AbstractGPUVector{SVector{N, Float32}},
    y::AbstractGPUVector{SVector{N, Float32}};
    groupsize::Int = 256
) where {N}

    #   Info
    space = SamplingSpace(core.shape, x, y)
    samples = get_num_samples(core.sampling, space)

    #   Allocate memory
    pv = get_power_vector(core, core.shape)
    offsets = get_offsets(core, core.shape)

    hist = KernelAbstractions.zeros(core.backend, Int32, get_histogram_size(core.shape))

    #   Call the kernel
    if core.sampling isa Full
        gpu_rng = KernelAbstractions.zeros(core.backend, Int32, 1)
        gpu_histogram!(core.backend, groupsize)(x, y, pv, offsets, core, space, Int32(samples), hist, gpu_rng, Int32(N); ndrange = samples)
    else
        rng = get_sample(core, core.sampling, space, samples)
        gpu_rng = KernelAbstractions.zeros(core.backend, SVector{2,Int32}, samples)
        copyto!(gpu_rng, rng)

        gpu_histogram!(core.backend, groupsize)(x, y, pv, offsets, core, space, Int32(samples), hist, gpu_rng, Int32(N); ndrange = samples)
    end

    KernelAbstractions.synchronize(core.backend)
    res =  hist |> Vector
    out = eachindex(res)

    return Counts(res, out)
end

##########################################################################################
#   Implementation: GPU Kernels
##########################################################################################
@kernel function gpu_histogram!(x, y, pv, offsets, core, space, samples, hist, rng, n)
    m = @index(Global)
    if m <= samples
        i = zero(Int32)
        j = zero(Int32)

        if core.sampling isa Full
            i, j = get_sample(core, core.sampling, space, rng, m)
        else
            i = rng[m][1]
            j = rng[m][2]
        end

        idx = gpu_compute_motif(core.shape.expr, x, y, i, j, pv, offsets, n)
        
        Atomix.@atomic hist[idx] += one(Int32)
    end
end

##########################################################################################
#   Implementation: distribution
##########################################################################################
distribution(
    x::AbstractGPUVector{SVector{N, Float32}},
    y::AbstractGPUVector{SVector{N, Float32}},
    shape::MicrostateShape;
    rate::Float32 = 0.05f0,
    sampling::SamplingMode = SRandom(rate),
    groupsize::Int = 256,
    backend = get_backend(x)
) where {N} = distribution(GPUCore(backend, shape, sampling), x, y; groupsize = groupsize)
#.........................................................................................
distribution(
    x::AbstractGPUVector{SVector{N, Float32}}, 
    y::AbstractGPUVector{SVector{N, Float32}},
    expr::RecurrenceExpression,
    n::Int;
    rate::Float32 = 0.05f0,
    sampling::SamplingMode = SRandom(rate),
    groupsize::Int = 256,
    backend = get_backend(x)
) where {N} = distribution(GPUCore(backend, Rect(expr, n), sampling), x, y; groupsize = groupsize)
#.........................................................................................
distribution(
    x::AbstractGPUVector{SVector{N, Float32}},
    y::AbstractGPUVector{SVector{N, Float32}},
    ε::Float32,
    n::Int;
    rate::Float32 = 0.05f0,
    sampling::SamplingMode = SRandom(rate),
    groupsize::Int = 256,
    backend = get_backend(x),
    metric::GPUMetric = GPUEuclidean()
) where {N} = distribution(x, y, Standard(ε; metric = metric), n; rate = rate, sampling = sampling, groupsize = groupsize, backend = backend)
#.........................................................................................
distribution(
    x::AbstractGPUVector{SVector{N, Float32}},
    shape::MicrostateShape;
    rate::Float32 = 0.05f0,
    sampling::SamplingMode = SRandom(rate),
    groupsize::Int = 256,
    backend = get_backend(x),
) where{N} = distribution(x, x, shape; rate = rate, sampling = sampling, groupsize = groupsize, backend = backend)
#.........................................................................................
distribution(
    x::AbstractGPUVector{SVector{N, Float32}},
    expr::RecurrenceExpression,
    n::Int;
    rate::Float32 = 0.05f0,
    sampling::SamplingMode = SRandom(rate),
    groupsize::Int = 256,
    backend = get_backend(x),
) where{N} = distribution(x, x, expr, n; rate = rate, sampling = sampling, groupsize = groupsize, backend = backend)
#.........................................................................................
distribution(
    x::AbstractGPUVector{SVector{N, Float32}},
    ε::Float32,
    n::Int;
    rate::Float32 = 0.05f0,
    sampling::SamplingMode = SRandom(rate),
    groupsize::Int = 256,
    backend = get_backend(x),
    metric::GPUMetric = GPUEuclidean()
) where {N} = distribution(x, x, ε, n; rate = rate, sampling = sampling, groupsize = groupsize, backend = backend, metric = metric)
#.........................................................................................
"""
    distribution(core::GPUCore, [x], [y]; kwargs...)

Compute an RMA distribution for the input data `[x]` and `[y]` using a GPU backend configuration
defined by `core`, which must be a [`GPUCore`](@ref).

The inputs `[x]` and `[y]` must be vectors of type `AbstractGPUVector`. This method supports
time-series analysis only.

!!! note
    The resulting distribution is copied from GPU memory back to the CPU.

### Arguments
- `core`: A [`GPUCore`](@ref) defining the [`MicrostateShape`](@ref),
  [`RecurrenceExpression`](@ref), and [`SamplingMode`](@ref) used in the computation.
- `[x]`: Input data provided as an `AbstractGPUVector`.
- `[y]`: Input data provided as an `AbstractGPUVector`.

### Keyword Arguments
- `groupsize`: Number of threads per GPU workgroup.

### Examples
```julia
using CUDA
gpudata = StateSpaceSet(Float32.(data)) |> CuVector
core = GPUCore(CUDABackend(), Rect(Standard(0.27f0; metric = GPUEuclidean()), 2), SRandom(0.05))
dist = distribution(core, gpudata, gpudata)
```

!!! warning
    Spatial data are not supported by [`GPUCore`](@ref).
"""
function distribution(
    core::GPUCore,
    x,
    y;
    groupsize::Int = 256,
)
    hist = histogram(core, x, y; groupsize = groupsize)
    return Probabilities(hist)
end

##########################################################################################