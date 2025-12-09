export Rect

##########################################################################################
#   MotifShape: Rect + Constructors and sub-types
##########################################################################################
abstract type Rect <: MotifShape end
#.........................................................................................
#   Based on time series: (CPU & GPU)
#.........................................................................................
struct Rect2{W, H, B, E <: RecurrenceExpression} <: Rect
    expr::E
end

Rect(expr::E; W = 2, H = 2, B = 2) where {E <: RecurrenceExpression} = Rect2{W,H,B,E}(expr)
Rect(expr::E, N; B = 2) where {E <: RecurrenceExpression} = Rect2{N,N,B,E}(expr)

#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
struct RectN{D, B, E <: RecurrenceExpression} <: Rect
    expr::E
    structure::NTuple{D, Int}
end

Rect(expr::E, structure::NTuple{D, Int}; B = 2) where {D, E <: RecurrenceExpression} = RectN{D, B, E}(expr, structure)

##########################################################################################
#   Implementations: SamplingSpace
##########################################################################################
#   Based on time series: (CPU & GPU)
#.........................................................................................
SamplingSpace(
    ::Rect2{W, H, B, E}, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{N, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{N, Float32}}}
) where {W, H, B, E<:RecurrenceExpression, N} = SSRect2(length(x) - W, length(y) - H)
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function SamplingSpace(
    shape::RectN{D, B, E}, 
    x::AbstractArray{<: Real}, 
    y::AbstractArray{<: Real}
) where {D, B, E<:RecurrenceExpression}

    dims_x = size(x)[2:end]
    dims_y = size(y)[2:end]

    dims = (dims_x..., dims_y...)

    @assert length(dims) == D "The motif shape and the input data dimension is not compatible."
    
    space = ntuple(i -> dims[i] - shape.structure[i], D)
    return SSRectN{D}(space)
end

##########################################################################################
#   Implementations: compute_motif (SRP)
##########################################################################################
@inline function compute_motif(
    shape::RectN,
    x::AbstractArray{<: Real},
    y::AbstractArray{<: Real},
    idx::Vector{Int},
    itr::Vector{Int},
    power_vector::SVector{N, Int}
) where {N}
    
    index = 0
    dim = ndims(x) - 1
    copy!(itr, idx)

    @inbounds @fastmath for p in power_vector

        i = ntuple(k -> itr[k], dim)
        j = ntuple(k -> itr[dim + k], length(shape.structure) - dim)

        index += recurrence(shape.expr, x, y, i, j) * p

        itr[1] += 1
        for k in 1:length(shape.structure) - 1
            if (itr[k] > idx[k] + (shape.structure[k] - 1))
                itr[k] = idx[k]
                itr[k + 1] += 1
            else
                break
            end
        end
    end

    return index + 1
end

##########################################################################################
#   Implementations: utils — histogram size, power vector, and offsets
##########################################################################################
@generated function get_histogram_size(::Rect2{W, H, B, E}) where {W, H, B, E}
    size = B^(W*H)
    return :( $size )
end

@generated function get_power_vector(::CPUCore, ::Rect2{W, H, B, E}) where {W, H, B, E}
    N = W * H
    expr = :(SVector{$N}( $([:(B^$i) for i in 0:(N-1)]... ) ))
    return expr
end

@generated function get_offsets(::CPUCore, ::Rect2{W, H, B, E}) where {W, H, B, E}
    N = W * H
    elems = [ :(SVector{2, Int}($w, $h)) for w in 0:(W - 1) for h in 0:(H - 1)]
    return :( SVector{$N, $(SVector{2, Int})}( $(elems...) ) )
end

@generated function get_power_vector(::GPUCore, ::Rect2{W, H, B, E}) where {W, H, B, E}
    N = W * H
    expr = :(SVector{$N}( $([:(Int32(B^$i)) for i in 0:(N-1)]... ) ))
    return expr
end

@generated function get_offsets(::GPUCore, ::Rect2{W, H, B, E}) where {W, H, B, E}
    N = W * H
    elems = [ :(SVector{2, Int32}($(Int32(w)), $(Int32(h)))) for w in 0:(W - 1) for h in 0:(H - 1)]
    return :( SVector{$N, $(SVector{2, Int32})}( $(elems...) ) )
end

function get_histogram_size(shape::RectN{D, B, E}) where {D, B, E}
    size = B^(prod(shape.structure))
    return size
end

function get_power_vector(::CPUCore, shape::RectN{D, B, E}) where {D, B, E}
    N = prod(shape.structure)
    return SVector{N}((B^i for i in 0:(N-1))...)
end