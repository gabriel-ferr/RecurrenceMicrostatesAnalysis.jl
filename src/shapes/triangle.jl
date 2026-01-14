export Triangle

##########################################################################################
#   MicrostateShape: Triangle + Constructors and sub-types
##########################################################################################
"""
    Triangle{N, B, E<:RecurrenceExpression} <: MicrostateShape

Triangle{N, B, E<:RecurrenceExpression} <: MicrostateShape

Define a triangular microstate shape, originally introduced by Hirata in 2021
[Hirata2021Triangle](@cite).

#   Constructor
```julia
Triangle(expr::E, N::Int; B::Int = 2)
```
where `expr` is the [`RecurrenceExpression`](@ref) used to evaluate recurrences and `N`
defines the size of the triangular microstate.

#   Example
```julia
n = 3
triangle = Triangle(expr, n)
```
The corresponding microstate structure is given by:
```math
\\begin{pmatrix}
R_{i,j} & R_{i,j + 1} & R_{i,j + 2} \\\\
 & R_{i + 1,j + 1} & R_{i + 1,j + 2} \\\\
 & & R_{i + 2,j + 2} \\\\
\\end{pmatrix}
```

!!! compat
    Triangular microstate shape is not compatible with spatial data.
"""
struct Triangle{N, B, E<:RecurrenceExpression} <: MicrostateShape
    expr::E
end

Triangle(expr::E, N::Int; B::Int = 2) where {E <: RecurrenceExpression} = Triangle{N,B,E}(expr)

##########################################################################################
#   Implementations: SamplingSpace
##########################################################################################
#   Based on time series: (CPU & GPU)
#.........................................................................................
SamplingSpace(
    ::Triangle{N, B, E}, 
    x::Union{StateSpaceSet, AbstractGPUVector{SVector{D, Float32}}}, 
    y::Union{StateSpaceSet, AbstractGPUVector{SVector{D, Float32}}}
) where {N, B, E<:RecurrenceExpression, D} = SSRect2(length(x) - N + 1, length(y) - N + 1)

##########################################################################################
#   Implementations: utils — histogram size, power vector, and offsets
##########################################################################################
@generated function get_histogram_size(::Triangle{N, B, E}) where {N, B, E}
    size = B^((N * (N + 1)) ÷ 2)
    return :( $size )
end

@generated function get_power_vector(::CPUCore, ::Triangle{N, B, E}) where {N, B, E}
    expr = :(SVector{$(N*(N + 1) ÷ 2)}( $([:(B^$( (((j - 1) * j) ÷ 2) + (i - 1) )) for j in 1:N for i in 1:j]... ) ))
    return expr
end

@generated function get_offsets(::CPUCore, ::Triangle{N, B, E}) where {N, B, E}
    elems = [ :(SVector{2, Int}($i, $j)) for j in 0:(N - 1) for i in 0:j]
    return :( SVector{$(N*(N + 1) ÷ 2), $(SVector{2, Int})}( $(elems...) ) )
end

@generated function get_power_vector(::GPUCore, ::Triangle{N, B, E}) where {N, B, E}
    expr = :(SVector{$(N*(N + 1) ÷ 2)}( $([:(Int32(B^$( (((j - 1) * j) ÷ 2) + (i - 1) ))) for j in 1:N for i in 1:j]... ) ))
    return expr
end

@generated function get_offsets(::GPUCore, ::Triangle{N, B, E}) where {N, B, E}
    elems = [ :(SVector{2, Int32}($(Int32(i)), $(Int32(j)))) for j in 0:(N - 1) for i in 0:j]
    return :( SVector{$(N*(N + 1) ÷ 2), $(SVector{2, Int32})}( $(elems...) ) )
end