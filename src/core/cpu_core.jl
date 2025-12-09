export CPUCore

##########################################################################################
#   RMACore: CPU
##########################################################################################
abstract type CPUCore{M<:MotifShape, S<:SamplingMode} end

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
function histogram(
    core::CPUCore,
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

    return reduce(+, fetch.(tasks))
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function histogram(
    core::CPUCore,
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

    return reduce(+, fetch.(tasks))
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

function distribution(
    core::CPUCore,
    x,
    y;
    threads = Threads.nthreads()
)
    hist = histogram(core, x, y; threads = threads)
    return Probabilities(hist)
end