export SamplingMode, SamplingSpace

##########################################################################################
#   Sampling Mode
##########################################################################################
"""
    SamplingMode

Abstract supertype defining how the initial position \$(i, j)\$ of each microstate is selected
during the construction of recurrence microstate distributions.

# Implementations
- [`SRandom`](@ref)
- [`Full`](@ref)
"""
abstract type SamplingMode end

##########################################################################################
#   Sampling Space
##########################################################################################
"""
    SamplingSpace

Define the range of valid indices used to sample the initial positions \$(i, j)\$ of microstates.
"""
abstract type SamplingSpace end
#.........................................................................................
#   Based on time series: RP & CRP (CPU)
#.........................................................................................
struct SSRect2 <: SamplingSpace
    W::Int
    H::Int
end
#.........................................................................................
SamplingSpace(
    shape::MicrostateShape, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{N, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{N, Float32}}}
) where {N} = throw("The sampling space is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
#.........................................................................................
#   Based on spatial data: SRP & CSRP (CPU)
#.........................................................................................
struct SSRectN{D} <: SamplingSpace
    space::NTuple{D, Int}
end
#.........................................................................................
SamplingSpace(
    shape::MicrostateShape, 
    x::AbstractArray{<: Real}, 
    y::AbstractArray{<: Real}
) = error("The sampling space is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")

##########################################################################################
#   Implementation: sampling
##########################################################################################
function get_sample(
    core::RMACore,
    mode::SamplingMode,
    space::SamplingSpace
)
    error("'get_sample' is not implemented for the set: \t\n Core: $(typeof(core)) \t\n Mode: $(typeof(mode)) \t\n Space: $(typeof(space))")
end

##########################################################################################
#   Utils: number of samples
##########################################################################################
function get_num_samples(
    mode::SamplingMode,
    space::SamplingSpace
)
    error("The number of samples is not implemented for a sampling space of type '$(typeof(space))' with sampling mode of type '$(typeof(mode))'")
end

##########################################################################################