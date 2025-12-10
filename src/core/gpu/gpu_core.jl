export GPUCore
using CUDA

##########################################################################################
#   RMACore: GPU
##########################################################################################
abstract type GPUCore{B, M<:MotifShape, S<:SamplingMode} end

struct StandardGPUCore{B, M<:MotifShape, S<:SamplingMode} <: GPUCore{B, M, S}
    backend::B
    shape::M
    sampling::S
end

GPUCore(backend::B, shape::M, sampling::S) where {B, M<:MotifShape, S<:SamplingMode} = StandardGPUCore(backend, shape, sampling)

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
function histogram(
    core::GPUCore,
    x::AbstractGPUVector{SVector{N, Float32}},
    y::AbstractGPUVector{SVector{N, Float32}};
    groupsize::Int = 256,
    seed::UInt32 = UInt32(time_ns() & 0xf12a57e8),
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
        gpu_rng = KernelAbstractions.zeros(core.backend, SVector{2,Int32}, 1)
        gpu_histogram!(core.backend, groupsize)(x, y, pv, offsets, core, space, Int32(samples), hist, gpu_rng, Int32(N); ndrange = samples)
    else
        rng = get_sample(core, core.sampling, space, samples)
        gpu_rng = KernelAbstractions.zeros(core.backend, SVector{2,Int32}, samples)
        copyto!(gpu_rng, rng)

        gpu_histogram!(core.backend, groupsize)(x, y, pv, offsets, core, space, Int32(samples), hist, gpu_rng, Int32(N); ndrange = samples)
    end

    KernelAbstractions.synchronize(core.backend)
    return hist |> Vector
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
            i, j = get_sample(core, core.sampling, space, nothing, m)
        else
            i = rng[m][1]
            j = rng[m][2]
        end

        idx = gpu_compute_motif(core.shape.expr, x, y, i, j, pv, offsets, n)

        if idx < 1 || idx > length(hist)
            @cuprintln("BAD IDX: ", idx)
        end

        
        Atomix.@atomic hist[idx] += one(Int32)
    end
end