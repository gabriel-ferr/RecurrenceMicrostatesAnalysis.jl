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

    K = 64  # quantos erros guardar
    badcount = KernelAbstractions.zeros(core.backend, Int32, 1)
    bad_m    = KernelAbstractions.zeros(core.backend, Int32, K)
    bad_i    = KernelAbstractions.zeros(core.backend, Int32, K)
    bad_j    = KernelAbstractions.zeros(core.backend, Int32, K)
    bad_idx  = KernelAbstractions.zeros(core.backend, Int32, K)
    bad_val  = KernelAbstractions.zeros(core.backend, Int32, K)  # valor de gpu_recurrence


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
        gpu_histogram!(core.backend, groupsize)(x, y, pv, offsets, core, space, Int32(samples), hist, gpu_rng, Int32(N),
        
        badcount, bad_m, bad_i, bad_j, bad_idx, bad_val
        
        ; ndrange = samples)
    else
        rng = get_sample(core, core.sampling, space, samples)
        gpu_rng = KernelAbstractions.zeros(core.backend, SVector{2,Int32}, samples)
        copyto!(gpu_rng, rng)

        gpu_histogram!(core.backend, groupsize)(x, y, pv, offsets, core, space, Int32(samples), hist, gpu_rng, Int32(N); ndrange = samples)
    end

    println("BEFORE SYNC")

    KernelAbstractions.synchronize(core.backend)

    println("AFTER SYNC")

    host_badcount = Vector{Int32}(undef,1); copyto!(host_badcount, badcount)
    if host_badcount[1] > 0
        c = host_badcount[1]
        host_m = Vector{Int32}(undef,K); copyto!(host_m, bad_m)
        host_i = Vector{Int32}(undef,K); copyto!(host_i, bad_i)
        host_j = Vector{Int32}(undef,K); copyto!(host_j, bad_j)
        host_ii = Vector{Int32}(undef,K); copyto!(host_ii, bad_ii)
        host_jj = Vector{Int32}(undef,K); copyto!(host_jj, bad_jj)
        host_idx = Vector{Int32}(undef,K); copyto!(host_idx, bad_idx)
        host_val = Vector{Int32}(undef,K); copyto!(host_val, bad_val)
        @show host_badcount[1]
        @show host_m[1:c], host_i[1:c], host_j[1:c], host_ii[1:c], host_jj[1:c], host_idx[1:c], host_val[1:c]
    end

    return hist |> Vector
end

##########################################################################################
#   Implementation: GPU Kernels
##########################################################################################
@kernel function gpu_histogram!(x, y, pv, offsets, core, space, samples, hist, rng, n,
     badcount, bad_m, bad_i, bad_j, bad_idx, bad_val)
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

        if idx < 1 || idx > length(hist) || i < 1 || i > space.W || j < 1 || j > space.H
            pos = Atomix.@atomic badcount[1] += 1
            if pos <= length(bad_m)
                bad_m[pos]  = m
                bad_i[pos]  = i
                bad_j[pos]  = j
                bad_idx[pos]= idx
                bad_val[pos]= Int32(val)  # val convertido
            end
        end
        
        Atomix.@atomic hist[idx] += one(Int32)
    end
end