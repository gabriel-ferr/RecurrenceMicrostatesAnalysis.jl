export Parameter, optimize

##########################################################################################
#   Parameter
##########################################################################################
"""
    Parameter

Abstract supertype for free parameters that can be optimized using RMA.

# Implementations
- [`Threshold`](@ref)
"""
abstract type Parameter end

##########################################################################################
#   Implementation: optimize
##########################################################################################
"""
    optimize(param::Parameter, qm::QuantificationMeasure, args...)

Optimize a free [`Parameter`](@ref) using the specified [`QuantificationMeasure`](@ref).

!!! warning
    The `optimize` function may compute multiple distributions and can be computationally expensive.
    Avoid calling it inside performance-critical loops.
"""
function optimize(param::Parameter, qm::QuantificationMeasure)
    error("The 'optimize' is not implemented to the parameter '$(typeof(param))' for the measure '$(typeof(qm))'.")
end

##########################################################################################