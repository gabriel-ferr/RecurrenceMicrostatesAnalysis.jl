export Diagonal

##########################################################################################
#   MotifShape: Diagonal + Constructors and sub-types
##########################################################################################
"""
    Diagonal{N, B, E<:RecurrenceExpression} <: MotifShape
"""
struct Diagonal{N, B, E<:RecurrenceExpression} <: MotifShape
    expr::E
end

Diagonal(expr::E, N::Int; B::Int = 2) where {E <: RecurrenceExpression} = Diagonal{N,B,E}(expr)

##########################################################################################
#   Implementations: SamplingSpace
##########################################################################################
#   Based on time series: (CPU & GPU)
#.........................................................................................
SamplingSpace(
    ::Diagonal{N, B, E}, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{D, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{D, Float32}}}
) where {N, B, E<:RecurrenceExpression, D} = SSRect2(length(x) - N, length(y) - N)

function SamplingSpace(
    ::Diagonal{N, B, E}, 
    x::AbstractArray{<: Real}, 
    y::AbstractArray{<: Real}
) where {N, B, E<:RecurrenceExpression}

    dims_x = size(x)[2:end]
    dims_y = size(y)[2:end]

    dims = (dims_x..., dims_y...)
    
    space = ntuple(i -> dims[i] - N, length(dims))
    return SSRectN{length(dims)}(space)
end

##########################################################################################
#   Implementations: compute_motif (SRP)
##########################################################################################
@inline function compute_motif(
    shape::Diagonal,
    x::AbstractArray{<: Real},
    y::AbstractArray{<: Real},
    idx::Vector{Int},
    itr::Vector{Int},
    power_vector::SVector{D, Int}
) where {D}
    
    index = 0
    dim_x = ndims(x) - 1
    dim_y = ndims(y) - 1

    copy!(itr, idx)

    @inbounds @fastmath for p in power_vector

        i = ntuple(k -> itr[k], dim_x)
        j = ntuple(k -> itr[dim_x + k], dim_y)

        index += recurrence(shape.expr, x, y, i, j) * p

        itr .+= 1
    end

    return index + 1
end

##########################################################################################
#   Implementations: utils — histogram size, power vector, and offsets
##########################################################################################
@generated function get_histogram_size(::Diagonal{N, B, E}) where {N, B, E}
    size = B^(N)
    return :( $size )
end

@generated function get_power_vector(::CPUCore, ::Diagonal{N, B, E}) where {N, B, E}
    expr = :(SVector{$N}( $([:(B^$i) for i in 0:(N-1)]... ) ))
    return expr
end

@generated function get_offsets(::CPUCore, ::Diagonal{N, B, E}) where {N, B, E}
    elems = [ :(SVector{2, Int}($w, $w)) for w in 0:(N - 1)]
    return :( SVector{$N, $(SVector{2, Int})}( $(elems...) ) )
end

@generated function get_power_vector(::GPUCore, ::Diagonal{N, B, E}) where {N, B, E}
    expr = :(SVector{$N}( $([:(Int32(B^$i)) for i in 0:(N-1)]... ) ))
    return expr
end

@generated function get_offsets(::GPUCore, ::Diagonal{N, B, E}) where {N, B, E}
    elems = [ :(SVector{2, Int32}($(Int32(w)), $(Int32(w)))) for w in 0:(W - 1)]
    return :( SVector{$N, $(SVector{2, Int32})}( $(elems...) ) )
end