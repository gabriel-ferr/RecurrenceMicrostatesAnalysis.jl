export QuantificationMeasure, measure

##########################################################################################
#   Type: Quantification Measure
##########################################################################################
"""
    QuantificationMeasure

Abstract supertype defining an RQA or RMA quantification measure.

All quantifiers implemented in the package subtype `QuantificationMeasure` and define their
computation via the [`measure`](@ref) function.

# Implementations
- [`Determinism`](@ref)
- [`Disorder`](@ref)
- [`Laminarity`](@ref)
- [`RecurrenceEntropy`](@ref)
- [`RecurrenceRate`](@ref)
"""
abstract type QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
"""
    measure(qm::QuantificationMeasure, [...])

Compute the quantification measure defined by the given [`QuantificationMeasure`](@ref) instance.

The accepted arguments (`[...]`) depend on the specific quantifier implementation.
"""
function measure(ms::QuantificationMeasure)
    error("There isn't a 'measure' implementation for '$(typeof(ms))'.")
end

##########################################################################################