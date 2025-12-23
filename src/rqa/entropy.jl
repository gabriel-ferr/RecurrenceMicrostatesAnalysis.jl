export RecurrenceEntropy

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
"""
    RecurrenceEntropy <: QuantificationMeasure

Defines the *Recurrence Microstates Entropy* (RME) quantification measure.

The computation of RME is performed using the [`measure`](@ref) function, for which two implementations are provided.

## Using a distribution
```julia
measure(::RecurrenceEntropy, dist::Probabilities)
```
### Input
- The `QuantificationMeasure`.
- `dist`: a distribution of recurrence microstates.

### Output
Returns a `Float64` corresponding to the RME computed using the Shannon entropy.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, 0.27, 3)
rme = measure(RecurrenceEntropy(), dist)
```

## Using a time series
```julia
measure(::RecurrenceEntropy, [x]; kwargs...)
```
### Input
- The `QuantificationMeasure`.
- `[x]`: time-series provided as an [`StateSpaceSet`](@ref).

### Output
Returns a `Float64` corresponding to the **maximum** RME based on the Shannon Entropy.

### Keyword arguments
- `n`: an `Integer` defining the microstates size. The default value is `3`.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
rme = measure(RecurrenceEntropy(), data; n = 4)
```
"""
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