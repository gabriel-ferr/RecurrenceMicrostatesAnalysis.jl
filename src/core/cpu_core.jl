export CPUCore

##########################################################################################
#   RMACore: CPU
##########################################################################################
struct CPUCore{M<:MotifShape, S<:SamplingMode} <: RMACore
    shape::M
    sampling::S
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
                idx = compute_motif(core.shape, x, y, i, j, pv, offsets)
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

    @assert ndims(x) - 1 + ndims(y) - 1 == length(core.shape.structure)
    @assert length(space.space) == length(core.shape.structure)

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