export MotifShape, compute_motif

##########################################################################################
#   MotifShape
##########################################################################################
abstract type MotifShape end

##########################################################################################
#   Index computation
##########################################################################################
#   Based on time series: (CPU & GPU)
#.........................................................................................
function compute_motif(
    shape::M,
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{DX, Float32}}},
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{DY, Float32}}},
    i::I,
    j::I,
    power_vector::SVector{N, I},
    offsets::SVector{N, I}
) where {M <: MotifShape, DX, DY, I <: Integer, N}
    throw("The index computation is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function compute_motif(
    shape::M,
    x::AbstractArray{<: Real},
    y::AbstractArray{<: Real},
    idx::Vector{Int},
    itr::Vector{Int},
    power_vector::SVector{N, Int}
) where {M <: MotifShape, N}
    throw("The index computation is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
##########################################################################################
#   Utils: number of recurrences and power vector.
##########################################################################################
@generated function get_histogram_size(
    shape::M
) where {M <: MotifShape}
    throw("The histogram size is not implemented for a motif shape of type '$(typeof(shape))'.")
end

@generated function get_power_vector(
    shape::M
) where {M <: MotifShape}
    throw("The power vector is not implemented for a motif shape of type '$(typeof(shape))'.")
end

@generated function get_offsets(
    core::C,
    shape::M
) where {C <: RMACore, M <: MotifShape}
    throw("Motif's offsets is not implemented for a motif shape of type '$(typeof(shape))' and a core of type '$(typeof(core))'.")
end