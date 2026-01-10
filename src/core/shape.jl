export MicrostateShape

##########################################################################################
#   MicrostateShape
##########################################################################################
"""
    MicrostateShape

Abstract supertype defining the basic structure and layout of a microstate (or motif).

A `MicrostateShape` specifies which relative positions are retrieved from the data to evaluate
recurrences, and how these binary recurrence values are interpreted and mapped to a decimal
representation for counting.

All subtypes of `MicrostateShape` must include a field `expr`, which defines the
[`RecurrenceExpression`](@ref) used to compute recurrences.

# Implementations
- [`Diagonal`](@ref)
- [`Rect`](@ref)
- [`Triangle`](@ref)
"""
abstract type MicrostateShape end

##########################################################################################
#   Index computation
##########################################################################################
function compute_motif()
    throw("The index computation is not implemented without arguments.")
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function compute_motif(
    shape::MicrostateShape,
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
function get_histogram_size(shape::MicrostateShape)
    error("The histogram size is not implemented for a motif shape of type '$(typeof(shape))'.")
end
#.........................................................................................
function get_power_vector(core::RMACore, shape::MicrostateShape)
    error("The power vector is not implemented for a motif shape of type '$(typeof(shape))' and a core of type '$(typeof(core))'.")
end
#.........................................................................................
function get_offsets(core::RMACore, shape::MicrostateShape)
    error("Motif's offsets is not implemented for a motif shape of type '$(typeof(shape))' and a core of type '$(typeof(core))'.")
end

##########################################################################################