export Laminarity

##########################################################################################
#   Quantification Measure: Laminarity
##########################################################################################
struct Laminarity <: QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
#       Using as input a RMA distribution.
#.........................................................................................
function measure(::Laminarity, dist::Probabilities)
    rr = measure(RecurrenceRate(), dist)
    return 1 - ((1/rr) * dist[3])
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::Laminarity, x::StateSpaceSet; threshold::Real = optimize(Threshold(), RecurrenceEntropy(), x, 3)[1])
    dist = distribution(x, Rect(Standard(threshold); W = 1, H = 3))
    measure(Laminarity(), dist)
end