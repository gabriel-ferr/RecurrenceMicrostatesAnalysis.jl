export RecurrenceExpression, recurrence

##########################################################################################
#   Recurrence Expression
##########################################################################################
"""
    RecurrenceExpression
"""
abstract type RecurrenceExpression end

##########################################################################################
#   Function: recurrence
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
function recurrence(
    expr::RecurrenceExpression,
    x::StateSpaceSet,
    y::StateSpaceSet,
    ::Int,
    ::Int,
)
    error("The recurrence computation is not implemented to a recurrence expression of type '$(typeof(expr))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end
#.........................................................................................
#   Based on spatial data: (CPU)
#.........................................................................................
function recurrence(
    expr::RecurrenceExpression,
    x::AbstractArray{<:Real},
    y::AbstractArray{<:Real},
    ::NTuple{N, Int}, 
    ::NTuple{M, Int},
) where {N, M} 
    error("The recurrence computation is not implemented to a recurrence expression of type '$(typeof(expr))' with input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end