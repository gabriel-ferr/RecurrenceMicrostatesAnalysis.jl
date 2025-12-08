export Corridor

##########################################################################################
#   RecurrenceExpression + Constructors
##########################################################################################
struct Corridor{F <: Real, M <: Metric} <: RecurrenceExpression
    ε_min::F
    ε_max::F
    metric::M
end

function Corridor(ε_min::Real, ε_max::Real; metric::Metric = DEFAULT_METRIC)
    return Corridor(ε_min, ε_max, metric)
end
##########################################################################################
#   Implementations
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
@inline function recurrence(
    expr::Corridor,
    x::StateSpaceSet,
    y::StateSpaceSet,
    i::Int,
    j::Int,
)
    distance = @inbounds evaluate(expr.metric, x.data[i], y.data[j])
    return UInt8(distance ≥ expr.ε_min && distance ≤ expr.ε_max)
end
#.........................................................................................
#   Based on time series: (GPU)
#.........................................................................................
@inline function gpu_recurrence(expr::Corridor, x, y, i, j, n)
    distance = gpu_evaluate(expr.metric, x, y, i, j, n)
    return UInt8(distance ≥ expr.ε_min && distance ≤ expr.ε_max)
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
@inline function recurrence(
    expr::Corridor,
    x::AbstractArray{<:Real},
    y::AbstractArray{<:Real},
    i::NTuple{N, Int}, 
    j::NTuple{M, Int},
) where {N, M} 
    distance = @inbounds evaluate(expr.metric, view(x, :, i...), view(y, :, j...))
    return UInt8(distance ≥ exp.ε_min && distance ≤ expr.ε_max)
end