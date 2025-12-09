export Triangle

##########################################################################################
#   MotifShape: Triangle + Constructors and sub-types
##########################################################################################
struct Triangle{N, B, E<:RecurrenceExpression} <: MotifShape
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
) where {N, B, E<:RecurrenceExpression, D} = SSRect2(length(x) - N, length(y) - N)

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