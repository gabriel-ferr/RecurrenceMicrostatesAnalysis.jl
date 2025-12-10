export Determinism

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
struct Determinism <: QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
#       Using as input a RMA distribution.
#.........................................................................................
function measure(::Determinism, dist::Probabilities)
    rr = measure(RecurrenceRate(), dist)
    return 1 - ((1/rr) * dist[3])
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::Determinism, x::StateSpaceSet; threshold::Real = 0.27, n::Integer = 3)
    dist = distribution(x, Diagonal(Standard(threshold), 3))
    measure(Determinism(), dist)
end