export QuantificationMeasure, measure

##########################################################################################
#   Type: Quantification Measure
##########################################################################################
"""
    QuantificationMeasure

Abstract supertype defining an RQA or RMA quantification measure.  
All quantifiers implemented in the package inherit from this type and define their computation through the [`measure`](@ref) function.
"""
abstract type QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
"""
    measure(qm::QuantificationMeasure, [...])

Compute the quantifier defined by the given [`QuantificationMeasure`](@ref) instance.  
Each implementation may accept different parameters (`[...]`), depending on the specific quantifier.
"""
function measure(ms::QuantificationMeasure)
    error("There isn't a 'measure' implementation for '$(typeof(ms))'.")
end