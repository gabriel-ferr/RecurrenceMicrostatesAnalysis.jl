export Operation, operate

##########################################################################################
#   Operations with microstates
##########################################################################################
"""
    Operation

Abstract supertype for operations that can be applied to recurrence microstates or to
recurrence microstate distributions.

#   Implementations:
- [`PermuteColumns`](@ref)
- [`PermuteRows`](@ref)
- [`Transpose`](@ref)
"""
abstract type Operation end

##########################################################################################
#   Function to operate
##########################################################################################
"""
    operate(op::Operation, [...])

Apply the operation defined by the given [`Operation`](@ref) instance.

The accepted arguments (`[...]`) depend on the specific operation implementation.
"""
function operate(op::Operation)
    T = typeof(op)
    msg = "`operate` not implemented without arguments for $T."
    throw(ArgumentError(msg))
end

##########################################################################################