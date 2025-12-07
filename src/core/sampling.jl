export SamplingMode, SamplingSpace

##########################################################################################
#   Sampling Mode
##########################################################################################
abstract type SamplingMode end

##########################################################################################
#   Sampling Space
##########################################################################################
abstract type SamplingSpace end
#.........................................................................................
#   Based on time series: RP & CRP
#.........................................................................................
struct SSRect2 <: SamplingSpace
    W::Int
    H::Int
end

SamplingSpace(
    shape::M, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{DX, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{DY, Float32}}}
) where {M <: MotifShape, DX, DY} = throw("The sampling space is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
#.........................................................................................
#   Based on spatial data: SRP & CSRP
#.........................................................................................
struct SSRectN{D} <: SamplingSpace
    space::NTuple{D, Int}
end

SamplingSpace(
    shape::M, 
    x::AbstractArray{<: Real}, 
    y::AbstractArray{<: Real}
) where {M <: MotifShape} = throw("The sampling space is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")

##########################################################################################
#   Implementation: sampling
##########################################################################################
function get_sample(
    core::C,
    mode::M,
    space::S,
    args...
) where {C <: RMACore, M <: SamplingMode, S <: SamplingSpace}
    throw("'get_sample' is not implemented for the set: \t\n Core: $(typeof(core)) \t\n Mode: $(typeof(mode)) \t\n Space: $(typeof(space))")
end

##########################################################################################
#   Utils: number of samples
##########################################################################################
function get_num_samples(
    mode::M,
    space::S
) where {M <: SamplingMode, S <: SamplingSpace}
    throw("The number of samples is not implemented for a sampling space of type '$(typeof(space))' with sampling mode of type '$(typeof(mode))'")
end