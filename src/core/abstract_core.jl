export RMACore, histogram, distribution

##########################################################################################
#   RMACore
##########################################################################################
"""
    RMACore

Abstract supertype that defines the execution pipeline of the package.

An instance of **RMACore** must be provided to the [`histogram`](@ref) function to determine how the
histogram computation is performed.

Concrete implementations of `RMACore` are [`CPUCore`](@ref) and [`GPUCore`](@ref), which target CPU
and GPU execution, respectively. Implementing custom subtypes of `RMACore` is **strongly discouraged**,
as doing so requires reimplementing several internal utilities for the package ecosystem to function
correctly.

#   Implementations
- [`CPUCore`](@ref)
- [`GPUCore`](@ref)
"""
abstract type RMACore end

##########################################################################################
#   Implementations: histogram & distribution
##########################################################################################
"""
    histogram(core::RMACore, [x], [y])

Compute the histogram of recurrence microstates for the input data `[x]` and `[y]` using the specified
backend `core`, which must be an [`RMACore`](@ref).

This function executes the full backend pipeline: sampling the recurrence space, constructing
microstates, and evaluating recurrences.

The result is returned as a [`Counts`](@ref) object, where each index corresponds to the decimal
representation of the associated microstate.
"""
function histogram(core::RMACore, x, y)
    core_type = typeof(core)
    x_type = typeof(x)
    y_type = typeof(y)

    msg = "`histogram` not implemented for $core_type and input data of types $x_type and $y_type."
    throw(ArgumentError(msg))
end
#.........................................................................................
"""
    distribution(core::RMACore, [x], [y])

Compute an RMA distribution from the recurrence structure constructed using the input data `[x]` and
`[y]`.

If `[x]` and `[y]` are identical, the result corresponds to a Recurrence Plot (RP); otherwise, it
corresponds to a Cross-Recurrence Plot (CRP).

The `core` argument must be a subtype of [`RMACore`](@ref) and defines how the computation is performed,
including the execution backend (CPU or GPU), the microstate shape, the recurrence expression, and the
sampling mode.

The result is returned as a [`Probabilities`](@ref) object, where each index corresponds to the decimal
representation of the associated microstate.

Internally, `distribution` calls [`histogram`](@ref) and converts the resulting counts into
probabilities.
"""
function distribution(core::RMACore, x, y)
    core_type = typeof(core)
    x_type = typeof(x)
    y_type = typeof(y)

    msg = "`distribution` not implemented for $core_type and input data of types $x_type and $y_type."
    throw(ArgumentError(msg))
end

##########################################################################################