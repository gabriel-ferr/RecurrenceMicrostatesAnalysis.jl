export Laminarity

##########################################################################################
#   Quantification Measure: Laminarity
##########################################################################################
"""
    Laminarity <: QuantificationMeasure

Defines the *Laminarity* (LAM) quantification measure.
The computation of LAM is performed using the [`measure`](@ref) function, for which two implementations are provided.

## Using a distribution
```julia
measure(::Laminarity, dist::Probabilities)
```

### Input
- The `QuantificationMeasure`.
- `dist`: a distribution of recurrence microstates. The distribution must be computed from square or line microstates with size 3.

### Output
Returns a `Float64` corresponding to the estimated laminarity.

### Examples
- Using **square** microstates:
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, 0.27, 3)
lam = measure(Laminarity(), dist)
```

- Using **line** microstates:
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, Rect(Standard(0.27); rows = 1, cols = 3))
lam = measure(Laminarity(), dist)
```

## Using a time series
```julia
measure(::Laminarity, [x]; kwargs...)
```

### Input
- The `QuantificationMeasure`.
- `[x]`: time-series data provided as a [`StateSpaceSet`](@ref).

### Output
Returns a `Float64` corresponding to the estimated laminarity.

### Keyword arguments
- `threshold`: threshold used to compute the RMA distribution. By default, this is the threshold that maximizes the RME.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
lam = measure(Laminarity(), data)
```

!!! note
    When a time series is provided as input, RecurrenceMicrostatesAnalysis.jl uses `line` microstates by default.
"""
struct Laminarity <: QuantificationMeasure end

##########################################################################################
#   Implementation: measure
##########################################################################################
#       Using as input a RMA distribution.
#.........................................................................................
function measure(::Laminarity, dist::Probabilities)
    if (length(dist) == 512)
        rr = measure(RecurrenceRate(), dist)

        values = zeros(Int, 64)
        v_idx = 1

        for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1, a5 in 0:1, a6 in 0:1
            I_1 = 2 + 8 * a1 + 16 * a2 + 32 * a3 + 64 * a4 + 128 * a5 + 256 * a6
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
        error("Laminarity must be computed using square or line motifs with n = 3.")
    end
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::Laminarity, x::StateSpaceSet; threshold::Real = optimize(Threshold(), RecurrenceEntropy(), x, 3)[1])
    dist = distribution(x, Rect(Standard(threshold); rows = 1, cols = 3))
    measure(Laminarity(), dist)
end