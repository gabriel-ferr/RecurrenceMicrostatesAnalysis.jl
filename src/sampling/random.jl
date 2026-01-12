export SRandom

##########################################################################################
#   Sampling Mode: SRandom
##########################################################################################
"""
    SRandom{F<:Real} <: SamplingMode

Sampling mode that randomly selects microstate positions \$(i, j)\$ within the
[`SamplingSpace`](@ref).

#   Constructors
```julia
SRandom(num_samples::Int)
SRandom(ratio::Union{Float32, Float64})
```

The sampling mode can be initialized either by specifying the exact number of microstates
to sample or by providing a ratio of the total number of possible microstates.

#   Examples
```julia
s = SRandom(1000)   # Specify the exact number of sampled microstates
s = SRandom(0.05)   # Specify a ratio of the total possible microstates
```
"""
struct SRandom{F <: Real} <: SamplingMode
    sampling_factor::F
end
#.........................................................................................
function SRandom(num_samples::Int)
    @assert num_samples ≥ 1 "The number of samples must be greater than 1."
    return SRandom{Int}(num_samples)
end 
#.........................................................................................
function SRandom(ratio::Union{Float32, Float64})
    @assert ratio > 0 "The sampling ratio must be greater than 0."
    @assert ratio ≤ 1.0 "The sampling ratio must be smaller than 1."

    return SRandom{typeof(ratio)}(ratio)
end

##########################################################################################
#   Implementation: sampling
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
function get_sample(::CPUCore, ::SRandom, space::SSRect2, rng, _)
    i = rand(rng, 1:space.W)
    j = rand(rng, 1:space.H)

    return i, j
end
#.........................................................................................
#   Based on time series: (GPU)
#.........................................................................................
function get_sample(::GPUCore, ::SRandom, space::SSRect2, samples)
    return [@SVector[Int32(rand(1:space.W)), Int32(rand(1:space.H))] for _ in 1:samples]
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function get_sample(::CPUCore, ::SRandom, space::SSRectN, idx::Vector{Int}, rng, _)
    for i in eachindex(idx)
        idx[i] = rand(rng, 1:space.space[i])
    end
end
##########################################################################################
#   Implementations: Utils
##########################################################################################
get_num_samples(mode::SRandom{<:Integer}, ::SSRect2) = mode.sampling_factor
get_num_samples(mode::SRandom{<:Real}, space::SSRect2) = ceil(Int, mode.sampling_factor * space.W * space.H)
#.........................................................................................
get_num_samples(mode::SRandom{<:Integer}, ::SSRectN) = mode.sampling_factor
get_num_samples(mode::SRandom{<:Real}, space::SSRectN) = ceil(Int, mode.sampling_factor * reduce(*, space.space))

##########################################################################################