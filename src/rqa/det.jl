export Determinism

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
"""
    Determinism <: QuantificationMeasure

Defines the *Determinism* (DET) quantification measure. 
The computation of DET is performed using the [`measure`](@ref) function, for which two implementations are provided.

## Using a distribution
```julia
measure(::Determinism, dist::Probabilities)
```

### Input
- The `QuantificationMeasure`.
- `dist`: a distribution of recurrence microstates. The distribution must be computed from **square** or **diagonal** microstates with size 3.

### Output
Returns a `Float64` corresponding to the estimated determinism.

### Examples
- Using **square** microstates:
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, 0.27, 3)
det = measure(Determinism(), dist)
```

- Using **diagonal** microstates:
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, Diagonal(Standard(0.27), 3))
det = measure(Determinism(), dist)
```

## Using a time series
```julia
measure(::Determinism, [x]; kwargs...)
```

### Input
- The `QuantificationMeasure`.
- `[x]`: time-series data provided as a [`StateSpaceSet`](@ref).

### Output
Returns a `Float64` corresponding to the estimated determinism.

### Keyword arguments
- `threshold`: threshold used to compute the RMA distribution. By default, this is the threshold that maximizes the RME.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
det = measure(Determinism(), data)
```

!!! note
    When a time series is provided as input, [`RecurrenceMicrostatesAnalysis`](@ref) uses [`Diagonal`](@ref) microstates by default.
"""
struct Determinism <: QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
#       Using as input a RMA distribution.
#.........................................................................................
function measure(::Determinism, dist::Probabilities)
    if (length(dist) == 512)
        rr = measure(RecurrenceRate(), dist)
        values = zeros(Int, 64)
        v_idx = 1

        for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1, a5 in 0:1, a6 in 0:1
            I_1 = 2 * a1 + 4 * a2 + 8 * a3 + 16 + 32 * a4 + 64 * a5 + 128 * a6
            values[v_idx] = I_1 + 1
            v_idx += 1
        end

        pl = 0.0
        for i in values
            pl += dist[i]
        end

        return 1 - ((1/rr) * pl)

    elseif (length(dist) == 8)
        rr = measure(RecurrenceRate(), dist)
        return 1 - ((1/rr) * dist[3])
    else
        error("Determinism must be computed using square or diagonal motifs with n = 3.")
    end
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::Determinism, x::StateSpaceSet; threshold::Real = optimize(Threshold(), RecurrenceEntropy(), x, 3)[1])
    dist = distribution(x, Diagonal(Standard(threshold), 3))
    measure(Determinism(), dist)
end