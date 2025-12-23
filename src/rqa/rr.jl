export RecurrenceRate

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
"""
    RecurrenceRate <: QuantificationMeasure

Defines the *Recurrence Rate* (RR) quantification measure. 
The computation of RR is performed using the [`measure`](@ref) function, for which two implementations are provided.

## Using a distribution
```julia
measure(::RecurrenceRate, dist::Probabilities)
```

### Input
- The `QuantificationMeasure`.
- `dist`: a distribution of recurrence microstates.

### Output
Returns a `Float64` corresponding to the estimated Recurrence Rate.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, 0.27, 3)
rr = measure(RecurrenceRate(), dist)
```

## Using a time series
```julia
measure(::RecurrenceRate, [x]; kwargs...)
```
### Input
- The `QuantificationMeasure`.
- `[x]`: time-series data provided as a [`StateSpaceSet`](@ref).

### Output
Returns a `Float64` corresponding to the estimated Recurrence Rate.

### Keyword arguments
- `n`: an `Integer` which defines the microstates size. The default value is `3`.
- `threshold`: threshold used to compute the RMA distribution. By default, this is the `threshold` that maximizes the RME.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
rme = measure(RecurrenceRate(), data; n = 4)
```
"""
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
        rr = count_ones(i - 1) / hv
        result += rr * dist[i]
    end

    return result
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::RecurrenceRate, x::StateSpaceSet; n::Integer = 3, threshold::Real = optimize(Threshold(), RecurrenceEntropy(), x, n)[1])
    dist = distribution(x, threshold, n)
    return measure(RecurrenceRate(), dist)
end