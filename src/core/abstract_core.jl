export RMACore, histogram, distribution

##########################################################################################
#   RMACore
##########################################################################################
"""
    RMACore
"""
abstract type RMACore end

##########################################################################################
#   Implementations: histogram & distribution
##########################################################################################
"""
    histogram(core::RMACore, [x], [y])

Compute the histogram of recurrence microstates using the input data `[x]` and `[y]` for a given `core`, which must be an [`RMACore`](@ref).  
This function implements the backend, executing the sampling process, constructing the microstates, and computing the recurrences.

The output is a [`Counts`](@ref) object, where each index corresponds to the decimal representation of the associated microstate.
"""
function histogram(core::RMACore, x, y)
    error("The RMA core of type '$(typeof(core))' is not implemented for input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end

"""
    distribution(core::RMACore, [x], [y])

Compute an RMA distribution from a recurrence structure constructed using the input data `[x]` and `[y]`.  
If `[x]` and `[y]` are identical, the result corresponds to a Recurrence Plot (RP); otherwise, it corresponds to a Cross-Recurrence Plot (CRP).

The `core` argument must be a structure inheriting from [`RMACore`](@ref) and defines how the computation is performed, including whether a CPU or GPU backend is used, the microstate shape, the recurrence expression, and the sampling mode.

The output of `distribution` is a [`Probabilities`](@ref) object, where each index corresponds to the decimal representation of the associated microstate.
"""
function distribution(core::RMACore, x, y)
    error("The RMA core of type '$(typeof(core))' is not implemented for input types: \n\t x: '$(typeof(x))'\n\t y: '$(typeof(y))')")
end

