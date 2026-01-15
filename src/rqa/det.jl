export Determinism

##########################################################################################
#   Quantification Measure: RecurrenceRate
##########################################################################################
"""
    Determinism <: QuantificationMeasure

Define the *Determinism* (DET) quantification measure.

DET can be computed either from a distribution of recurrence microstates or directly from
time-series data. In both cases, the computation is performed via the [`measure`](@ref)
function.

#   Using a distribution
```julia
measure(::Determinism, dist::Probabilities)
```

##  Arguments
- `dist`: A distribution of recurrence microstates. The distribution must be computed from
    **square** or **diagonal** microstates of size 3.

##  Returns
A `Float64` corresponding to the estimated determinism.

##  Examples
### Using square microstates
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, 0.27, 3)
det = measure(Determinism(), dist)
```

### Using diagonal microstates
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, Diagonal(Standard(0.27), 3))
det = measure(Determinism(), dist)
```

#   Using a time series
```julia
measure(::Determinism, [x]; kwargs...)
```

##  Arguments
- `[x]`: Time-series data provided as a [`StateSpaceSet`](@ref).

##  Returns
A `Float64` corresponding to the estimated determinism.

##  Keyword Arguments
- `threshold`: Threshold used to compute the RMA distribution. By default, this is chosen as
    the threshold that maximizes the recurrence microstate entropy (RME).

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
det = measure(Determinism(), data)
```

!!! note
    When time-series data are provided directly, RecurrenceMicrostatesAnalysis.jl uses [`Diagonal`](@ref) microstates by default.
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
        msg = "Determinism must be computed using square or diagonal microstates with n = 3."
        throw(ArgumentError(msg))
    end
end
#.........................................................................................
#       Using as input a time series
#.........................................................................................
function measure(::Determinism, x::StateSpaceSet; threshold::Real = optimize(Threshold(), RecurrenceEntropy(), x, 3)[1])
    dist = distribution(x, Diagonal(Standard(threshold), 3))
    measure(Determinism(), dist)
end

##########################################################################################