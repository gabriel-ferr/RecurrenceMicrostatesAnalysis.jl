export RecurrenceRate

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
struct RecurrenceRate <: QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
#       Using as input a RMA distribution.
#.........................................................................................
function measure(::RecurrenceRate, dist::Probabilities)
    result = 0.0
    hv = Int(log2(length(dist)))

    for i in eachindex(dist)
        rr = sum(digits(i - 1, base = 2)) / hv
        result += rr * dist[i]
    end

    return result
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::RecurrenceRate, x::StateSpaceSet; threshold::Real = 0.27, n::Integer = 3)
    dist = distribution(x, threshold, n)
    return measure(RecurrenceRate(), dist)
end