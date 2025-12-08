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
    shape::MotifShape,
    x::StateSpaceSet,
    y::StateSpaceSet,
    ::Int,
    ::Int,
    power_vector::SVector{N, Int},
    offsets::SVector{N, Int}
) where {N}
    throw("The index computation is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function compute_motif(
    shape::MotifShape,
    x::AbstractArray{<: Real},
    y::AbstractArray{<: Real},
    ::Vector{Int},
    ::Vector{Int},
    ::SVector{N, Int}
) where {N}
    error("The index computation is not implemented for a motif shape of type '$(typeof(shape))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
##########################################################################################
#   Utils: number of recurrences and power vector.
##########################################################################################
function get_histogram_size(shape::MotifShape)
    error("The histogram size is not implemented for a motif shape of type '$(typeof(shape))'.")
end

function get_power_vector(core::RMACore, shape::MotifShape)
    error("The power vector is not implemented for a motif shape of type '$(typeof(shape))' and a core of type '$(typeof(core))'.")
end

function get_offsets(core::RMACore, shape::MotifShape)
    error("Motif's offsets is not implemented for a motif shape of type '$(typeof(shape))' and a core of type '$(typeof(core))'.")
end