export Standard

##########################################################################################
#   RecurrenceExpression + Constructors
##########################################################################################
struct Standard{F <: Real, M <: Metric} <: RecurrenceExpression
    ε::F
    metric::M
end

function Standard(ε::Real; metric::Metric = DEFAULT_METRIC)
     return Standard(ε, metric)
end

##########################################################################################
#   Implementations
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
@inline function recurrence(
    expr::Standard,
    x::StateSpaceSet,
    y::StateSpaceSet,
    i::Int,
    j::Int,
)

    distance = @inbounds evaluate(expr.metric, x.data[i], y.data[j])
    return UInt8(distance ≤ expr.ε)
end
#.........................................................................................
#   Based on time series: (GPU)
#.........................................................................................
@inline function recurrence(
    expr::Standard,
    x::AbstractGPUVector{SVector{DX, Float32}},
    y::AbstractGPUVector{SVector{DY, Float32}},
    i::Int32,
    j::Int32,
) where {DX, DY}
    
    @inbounds xi = x[i]
    @inbounds yj = y[j]

    distance = @inbounds evaluate(expr.metric, xi, yj)
    return UInt8(distance ≤ expr.ε)
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
@inline function recurrence(
    expr::Standard,
    x::AbstractArray{<:Real},
    y::AbstractArray{<:Real},
    i::NTuple{N, Int}, 
    j::NTuple{M, Int},
) where {N, M} 
    
    distance = @inbounds evaluate(expr.metric, view(x, :, i...), view(y, :, j...))
    return UInt8(distance ≤ expr.ε)
end