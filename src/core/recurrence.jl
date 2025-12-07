export RecurrenceExpression, recurrence

##########################################################################################
#   Recurrence Expression
##########################################################################################
abstract type RecurrenceExpression end

##########################################################################################
#   Function: recurrence
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
function recurrence(
    expr::E,
    x::StateSpaceSet,
    y::StateSpaceSet,
    i::Int,
    j::Int,
) where {E <: RecurrenceExpression}
    throw("The recurrence computation is not implemented to a recurrence expression of type '$(typeof(expr))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
#.........................................................................................
#   Based on time series: (GPU)
#.........................................................................................
function recurrence(
    expr::E,
    x::AbstractGPUVector{SVector{DX, Float32}},
    y::AbstractGPUVector{SVector{DY, Float32}},
    i::Int32,
    j::Int32,
) where {E <: RecurrenceExpression, DX, DY}
    throw("The recurrence computation is not implemented to a recurrence expression of type '$(typeof(expr))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
#.........................................................................................
#   Based on spatial data: (CPU only)
#.........................................................................................
function recurrence(
    expr::E,
    x::AbstractArray{<:Real},
    y::AbstractArray{<:Real},
    i::NTuple{N, Int}, 
    j::NTuple{M, Int},
) where {E <: RecurrenceExpression, N, M} 
    throw("The recurrence computation is not implemented to a recurrence expression of type '$(typeof(expr))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end