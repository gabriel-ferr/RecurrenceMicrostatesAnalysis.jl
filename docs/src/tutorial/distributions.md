#   Distributions

This section introduces the computation of Recurrence Microstates Analysis (RMA) distributions
using **RecurrenceMicrostatesAnalysis.jl**.

We begin with a [Quick start with RecurrenceMicrostatesAnalysis.jl](@ref), which demonstrates a
simple application example. Next, we present [A brief review](@ref) of Recurrence Plots (RP)
and RMA. Finally, we explain the [`distribution`](@ref) function in
[Computing RMA distributions](@ref), including the role of [Histograms](@ref).

##  Quick start with RecurrenceMicrostatesAnalysis.jl

This section presents concise examples illustrating how to use the package. RMA distributions
are computed using the [`distribution`](@ref) function, which returns a [`Probabilities`](@ref)
structure containing the microstate distribution.

We start with a simple example based on a uniform random process. First, we generate the data
and convert it into a [`StateSpaceSet`](@ref):
```@example quick_example
using Distributions, RecurrenceMicrostatesAnalysis
data = rand(Uniform(0, 1), 10_000);
ssset = StateSpaceSet(data)
```

Next, we compute the RMA distribution. This requires specifying the recurrence threshold
$\varepsilon$ and the microstate size $N$. These parameters are discussed in more detail in [A brief review](@ref) and [Optimizing a parameter](@ref).

```@example quick_example
ε = 0.27
N = 2
dist = distribution(ssset, ε, N)
```

As another example, we use [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/)
to generate data from the Hénon map, following the example presented in its documentation:
```@example quick_example
using DynamicalSystems

function henon_rule(u, p, n) # here `n` is "time", but we don't use it.
    x, y = u # system state
    a, b = p # system parameters
    xn = 1.0 - a*x^2 + y
    yn = b*x
    return SVector(xn, yn)
end

u0 = [0.2, 0.3]
p0 = [1.4, 0.3]

henon = DeterministicIteratedMap(henon_rule, u0, p0)

total_time = 10_000
X, t = trajectory(henon, total_time)
X
```

Finally, we compute the RMA distribution of the trajectory `X`. Here, the threshold is selected
using [`optimize`](@ref) by maximizing the recurrence microstate entropy:
```@example quick_example
ε, S = optimize(Threshold(), RecurrenceEntropy(), X, N)
```

```@example quick_example
dist = distribution(X, ε, N)
```

##  A brief review

Recurrence Plots (RPs) were introduced in 1987 by Eckmann et al. [Eckmann1987RP](@cite) as a method for analyzing dynamical systems through recurrence
properties.

Consider a time series $\vec{x}_i \in \mathbb{R}^d$, $i \in \{1, 2, \dots, K\}$, where $K$ is the length of the time series and $d$ is the dimension of the phase space.
The recurrence plot is defined by the recurrence matrix
```math
R_{i,j} = \Theta(\varepsilon - \|\vec x_i - \vec x_j\|),
```
where $\Theta(\cdot)$ denotes the Heaviside step function and $\varepsilon$ is the recurrence
threshold.

The following figure shows examples of recurrence plots for different systems:
(a) white noise;
(b) a superposition of harmonic oscillators;
(c) a logistic map with a linear trend;
(d) Brownian motion.

![Image of four RPs with their timeseries](../assets/rps.png)

A recurrence microstate is a local structure extracted from an RP. For a given microstate
shape and size, the set of possible microstates is finite. For example, square microstates
with size $N = 2$ yield $16$ distinct configurations.

![Image of the 16 squared microstates to N = 2](../assets/microstates.png)

Recurrence Microstates Analysis (RMA) uses the probability distribution of these microstates
as a source of information for characterizing dynamical systems.

##  Computing RMA distributions

The computation of RMA distributions is the core functionality of
RecurrenceMicrostatesAnalysis.jl. All other tools in the package rely on these
distributions as their primary source of information.

RMA distributions are computed using the [`distribution`](@ref) function:
```@docs
distribution
```

A commonly used convenience interface is:
```julia
distribution([x], ε::Float, n::Int; kwargs...)
```
This method automatically selects a [`CPUCore`](@ref) when `x` is a [`StateSpaceSet`](@ref)
and a [`GPUCore`](@ref) when `x` is an `AbstractGPUVector`. By default, square microstates of size `n` are used.

Additional keyword arguments include:
- `rate::Float64`: Sampling rate (default: `0.05`).
- `sampling::SamplingMode`: Sampling mode (default: [`SRandom`](@ref)).
- `metric::Metric`: Distance metric from [Distances.jl](https://github.com/JuliaStats/Distances.jl). When using a [`GPUCore`](@ref), a [`GPUMetric`](@ref) must be provided.

!!! warning
    GPU backends require inputs of type `Float32`. `Float64` inputs are not supported on GPU.

Alternatively, a [`RecurrenceExpression`](@ref) can be specified directly:
```julia
distribution([x], expr::RecurrenceExpression, n::Int; kwargs...)
```
**Example:**
```@example quick_example
expr = Corridor(0.05, 0.27)
dist = distribution(ssset, expr, 2)
```

If a custom [`MicrostateShape`](@ref) is required, the call simplifies to:
```julia
distribution([x], shape::MicrostateShape; kwargs...)
```
**Example:**
```@example quick_example
shape = Triangle(Standard(0.27), 3)
dist = distribution(ssset, shape)
```

---
## Cross-recurrence plots
RMA distributions can also be computed from Cross-Recurrence Plots (CRPs) by providing two time series:
```julia
distribution([x], [y], expr::RecurrenceExpression, n::Int; kwargs...)
distribution([x], [y], shape::MicrostateShape; kwargs...)
```

**Example:**
```@example quick_example
data_1 = StateSpaceSet(rand(Uniform(0, 1), 1000))
data_2 = StateSpaceSet(rand(Uniform(0, 1), 2000))
dist = distribution(data_1, data_2, 0.27, 2)
```

!!! danger
    The inputs `x` and `y` must have the same phase-space dimensionality. The following example
    is invalid and will raise an exception:
    ```julia
    data_1 = StateSpaceSet(rand(Uniform(0, 1), (1000, 2)))
    data_2 = StateSpaceSet(rand(Uniform(0, 1), (2000, 3)))
    dist = distribution(data_1, data_2, 0.27, 2)
    ```

---
## Spatial data

The package also provides experimental support for spatial data, following *"Generalised Recurrence Plot Analysis for Spatial Data"* [Marwan2007Spatial](@cite).
In this context, input data are provided as `AbstractArray`s:
```math
    \vec{x}_{\vec i} \in \mathbb{R}^m,\quad \vec{i} \in \mathbb{Z}^d
```

For example:
```@example quick_example
spatialdata = rand(Uniform(0, 1), (2, 50, 50))
```
Due to the high dimensionality of spatial recurrence plots, direct visualization is often
infeasible. RMA distributions provide a compact alternative by sampling microstates directly
from the data.

**Examples:**
- Full $2 \times d$ microstates:
```@example quick_example
distribution(spatialdata, Rect(Standard(0.27), (2, 2, 2, 2)))
```
```@example quick_example
spatialdata_1 = rand(Uniform(0, 1), (2, 50, 50))
spatialdata_2 = rand(Uniform(0, 1), (2, 25, 25))
distribution(spatialdata_1, spatialdata_2, Rect(Standard(0.27), (2, 2, 2, 2)))
```

- Projected microstates:
```@example quick_example
distribution(spatialdata, Rect(Standard(0.27), (2, 1, 2, 1)))
```
```@example quick_example
distribution(spatialdata_1, spatialdata_2, Rect(Standard(0.27), (2, 1, 2, 1)))
```

##  Histograms
The [`histogram`](@ref) function counts the occurrences of each microstate identified during
sampling. It is called internally by [`distribution`](@ref) , which converts the resulting [`Counts`](@ref) into [`Probabilities`](@ref).

```@docs
histogram
```