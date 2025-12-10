export RecurrenceEntropy

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
struct RecurrenceEntropy <: QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
#       Using as input a RMA distribution.
#.........................................................................................
function measure(::RecurrenceEntropy, dist::Probabilities)
    return entropy(dist)
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::RecurrenceEntropy, x::StateSpaceSet; threshold::Real = 0.27, n::Integer = 3)
    dist = distribution(x, threshold, n)
    return measure(RecurrenceEntropy(), dist)
end