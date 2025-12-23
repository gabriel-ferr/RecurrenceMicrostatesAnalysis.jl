export CPUCore, StandardCPUCore

##########################################################################################
#   RMACore: CPU
##########################################################################################
"""
    CPUCore{M<:MotifShape, S<:SamplingMode}
"""
abstract type CPUCore{M<:MotifShape, S<:SamplingMode} end

"""
    StandardCPUCore{M<:MotifShape, S<:SamplingMode} <: CPUCore{M, S}
"""
struct StandardCPUCore{M<:MotifShape, S<:SamplingMode} <: CPUCore{M, S}
    shape::M
    sampling::S
end

CPUCore(shape::M, sampling::S) where {M<:MotifShape, S<:SamplingMode} = StandardCPUCore(shape, sampling)

##########################################################################################
#   Implementation: compute_motif
##########################################################################################
@inline function compute_motif(
    expr::E,
    x::StateSpaceSet,
    y::StateSpaceSet,
    i::Int,
    j::Int,
    power_vector::SVector{D, Int},
    offsets::SVector{D, SVector{2, Int}}
) where {E<:RecurrenceExpression, D}
    @inbounds begin
        index = 0

        for m in eachindex(power_vector)
            dw, dh = offsets[m]
            @fastmath index += recurrence(expr, x, y, i + dw, j + dh) * power_vector[m]
        end

        return index + 1
    end
end

##########################################################################################
#   Implementation: histogram
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
"""
    histogram(core::StandardCPUCore, [x], [y]; kwargs...)

Count the microstates of an "abstract" RP constructed using `[x]` and `[y]`.
If `[x]` and `[y]` are identical, the result corresponds to a Recurrence Plot (RP); otherwise, it corresponds to a Cross-Recurrence Plot (CRP).
The output is a histogram of recurrence microstates for the given input data as a [`Counts`](@ref) structure.

This method implements the CPU backend, based on a [`CPUCore`](@ref), specifically a [`StandardCPUCore`](@ref).

### Input
- `core`: A [`StandardCPUCore`](@ref), which defines how the backend computation is performed.
- `[x]`: Input data, provided as a [`StateSpaceSet`](@ref) or an `AbstractArray`.
- `[y]`: Input data, provided as a [`StateSpaceSet`](@ref) or an `AbstractArray`.

!!! note
    [`StateSpaceSet`](@ref) and `AbstractArray` use different backends and therefore different internal histogram implementations.  
    However, both functions share the same method signature, differing only in the input data format.

### Keyword Arguments
- `threads`: Number of threads used to compute the histogram. By default, this is set to `Threads.nthreads()`, which can be specified at Julia startup using -- using `--threads N` or via the environment variable `JULIA_NUM_THREADS`.

### Examples
- Time series:
```julia
core = CPUCore(Rect(Standard(0.27), 2), SRandom(0.05))
dist = histogram(core, ssset, ssset)
```

- Spatial data:
```julia
spatialdata = rand(Uniform(0, 1), (3, 50, 50))
core = CPUCore(Rect(Standard(0.5), (2, 2, 1, 1)), SRandom(0.05))
dist = histogram(core, spatialdata, spatialdata)
```
"""
function histogram(
    core::StandardCPUCore,
    x::StateSpaceSet,
    y::StateSpaceSet;
    threads = Threads.nthreads()
)
    #   Info
    space = SamplingSpace(core.shape, x, y)
    samples = get_num_samples(core.sampling, space)

    #   Allocate memory
    pv = get_power_vector(core, core.shape)
    offsets = get_offsets(core, core.shape)

    #   Compute the histogram
    chunk = ceil(Int, samples / threads)
    tasks = Vector{Task}(undef, threads)

    for t in 1:threads
        tasks[t] = Threads.@spawn begin
            local_hist = zeros(Int, get_histogram_size(core.shape))
            local_rng = TaskLocalRNG()

            start = (t - 1) * chunk + 1
            stop  = min(t * chunk, samples)

            for m in start:stop
                i, j = get_sample(core, core.sampling, space, local_rng, m)
                idx = compute_motif(core.shape.expr, x, y, i, j, pv, offsets)
                @inbounds local_hist[idx] += 1
            end

            return local_hist
        end
    end

    res = reduce(+, fetch.(tasks))
    out = eachindex(res)

    return Counts(res, out)
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function histogram(
    core::StandardCPUCore,
    x::AbstractArray{<: Real},
    y::AbstractArray{<: Real};
    threads = Threads.nthreads()
)
    #   Info
    space = SamplingSpace(core.shape, x, y)
    samples = get_num_samples(core.sampling, space)
    dim_x = ndims(x) - 1
    dim_y = ndims(y) - 1

    #   Allocate memory
    pv = get_power_vector(core, core.shape)

    #   Compute the histogram
    chunk = ceil(Int, samples / threads)
    tasks = Vector{Task}(undef, threads)

    for t in 1:threads
        tasks[t] = Threads.@spawn begin
            local_hist = zeros(Int, get_histogram_size(core.shape))
            local_rng = TaskLocalRNG()

            start = (t - 1) * chunk + 1
            stop  = min(t * chunk, samples)

            idx = zeros(Int, dim_x + dim_y)
            itr = zeros(Int, dim_x + dim_y)

            for m in start:stop
                get_sample(core, core.sampling, space, idx, local_rng, m)
                i = compute_motif(core.shape, x, y, idx, itr, pv)
                @inbounds local_hist[i] += 1
            end

            return local_hist
        end
    end

    res =  reduce(+, fetch.(tasks))
    out = eachindex(res)

    return Counts(res, out)
end

##########################################################################################
#   Implementation: distribution
##########################################################################################

distribution(
    x::StateSpaceSet,
    y::StateSpaceSet,
    shape::MotifShape;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads()
) = distribution(CPUCore(shape, sampling), x, y; threads = threads)

distribution(
    x::StateSpaceSet, 
    y::StateSpaceSet,
    expr::RecurrenceExpression,
    n::Int;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads()
) = distribution(CPUCore(Rect(expr, n), sampling), x, y; threads = threads)

distribution(
    x::StateSpaceSet,
    y::StateSpaceSet,
    ε::Float64,
    n::Int;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads(),
    metric::Metric = DEFAULT_METRIC
) = distribution(x, y, Standard(ε; metric = metric), n; rate = rate, sampling = sampling, threads = threads)

distribution(
    x::StateSpaceSet,
    shape::MotifShape;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads()
) = distribution(x, x, shape; rate = rate, sampling = sampling, threads = threads)

distribution(
    x::StateSpaceSet,
    expr::RecurrenceExpression,
    n::Int;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads()
) = distribution(x, x, expr, n; rate = rate, sampling = sampling, threads = threads)

distribution(
    x::StateSpaceSet,
    ε::Float64,
    n::Int;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads(),
    metric::Metric = DEFAULT_METRIC
) = distribution(x, x, ε, n; rate = rate, sampling = sampling, threads = threads, metric = metric)

distribution(
    x::AbstractArray{<:Real},
    y::AbstractArray{<:Real},
    shape::MotifShape;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads()
) = distribution(CPUCore(shape, sampling), x, y; threads = threads)

distribution(
    x::AbstractArray{<:Real},
    shape::MotifShape;
    rate::Float64 = 0.05,
    sampling::SamplingMode = SRandom(rate),
    threads::Int = Threads.nthreads()
) = distribution(x, x, shape; rate = rate, sampling = sampling, threads = threads)

"""
    distribution(core::CPUCore, [x], [y]; kwargs...)

Computes an RMA distribution using `[x]` and `[y]` as input and a CPU backend configuration specified by `core`, which must be a [`CPUCore`](@ref). 
The inputs `[x]` and `[y]` must be provided as [`StateSpaceSet`](@ref) objects for time-series data, or as `AbstractArray`s for spatial data.

### Input
- `core`: A [`CPUCore`](@ref) defining the configuration of the [`MotifShape`](@ref), [`RecurrenceExpression`](@ref), and [`SamplingMode`](@ref).
- `[x]`: Input data, given as a [`StateSpaceSet`](@ref) or an `AbstractArray`.
- `[y]`: Input data, given as a [`StateSpaceSet`](@ref) or an `AbstractArray`.

### Keyword Arguments
- `threads`: Number of threads which used to compute the distribution. By default, this is set to `Threads.nthreads()`, which can be specified at Julia startup using -- using `--threads N` or via the environment variable `JULIA_NUM_THREADS`.

### Examples
- Time series:
```julia
core = CPUCore(Rect(Standard(0.27), 2), SRandom(0.05))
dist = distribution(core, ssset, ssset)
```

- Spatial data:
```julia
spatialdata = rand(Uniform(0, 1), (3, 50, 50))
core = CPUCore(Rect(Standard(0.5), (2, 2, 1, 1)), SRandom(0.05))
dist = distribution(core, spatialdata, spatialdata)
```
"""
function distribution(
    core::CPUCore,
    x,
    y;
    threads = Threads.nthreads()
)
    hist = histogram(core, x, y; threads = threads)
    return Probabilities(hist)
end