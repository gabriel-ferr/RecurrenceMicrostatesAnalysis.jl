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
    T = typeof(param)
    msg = "`optimize` not implemented without arguments for $T using the quantification measure $(typeof(qm))."
    throw(ArgumentError(msg))
end

##########################################################################################