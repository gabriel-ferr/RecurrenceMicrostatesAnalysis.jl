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
function measure(::RecurrenceEntropy, x::StateSpaceSet; n::Integer = 3)
    return optimize(Threshold(), RecurrenceEntropy(), x, n)[2]
end