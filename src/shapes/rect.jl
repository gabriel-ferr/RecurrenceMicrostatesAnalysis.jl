export Rect

##########################################################################################
#   MicrostateShape: Rect + Constructors and sub-types
##########################################################################################
"""
    Rect <: MicrostateShape

Define a rectangular microstate shape.

`Rect` can represent either a two-dimensional microstate (identified as `Rect2`, used for
Recurrence Plots and Cross-Recurrence Plots) or an N-dimensional microstate (identified as
`RectN`, used for spatial data).

#   Rect2 (time-series data)
A 2D rectangular microstate can be initialized using either of the following constructors:
```julia
Rect(expr::E; rows = 2, cols = 2, B = 2)
Rect(rows::Int, cols::Int; expr = Standard(0.27), B = 2)
```

Here, `rows` and `columns` define the rectangle dimensions, `expr` is the [`RecurrenceExpression`](@ref)
used to evaluate recurrences, and `B` is the base used to encode microstate elements (typically `2`, representing recurrence or non-recurrence).

Rectangular microstates can be specialized to define common patterns such as lines,
columns, and squares:
```julia
line = Rect(expr; rows = n, cols = 1)
column = Rect(expr; rows = 1, cols = n)
square = Rect(expr; rows = n, cols = n)
```

Since square microstates are frequently used, a convenience constructor is also provided:
```julia
Rect(expr::E, N; B = 2)
```
where `N` defines the size of the square microstate. For example:
```julia
square = Rect(expr, n)
```

#   RectN (spatial data)
For N-dimensional structures, typically used with spatial data, the RectN variant can be
initialized as:
```julia
Rect(expr::E, structure::NTuple{D, Int}; B = 2)
```
Here, `structure` defines the size of the microstate along each dimension. For example:
```julia
nrect = Rect(expr, (2, 1, 2, 1))
```
This form is suitable for N-dimensional spatial data, such as images or volumetric datasets.
"""
abstract type Rect <: MicrostateShape end
#.........................................................................................
#   Based on time series: (CPU & GPU)
#.........................................................................................
struct Rect2{W, H, B, E <: RecurrenceExpression} <: Rect
    expr::E
end

Rect(expr::E; rows = 2, cols = 2, B = 2) where {E <: RecurrenceExpression} = Rect2{rows,cols,B,E}(expr)
Rect(expr::E, N; B = 2) where {E <: RecurrenceExpression} = Rect2{N,N,B,E}(expr)
Rect(rows::Int, cols::Int; expr::E = Standard(0.27), B = 2) where {E <: RecurrenceExpression} = Rect2{rows, cols, B, Standard}(expr)

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
) where {W, H, B, E<:RecurrenceExpression, N} = SSRect2(length(x) - W + 1, length(y) - H + 1)
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
    
    space = ntuple(i -> dims[i] - shape.structure[i] + 1, D)
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