#   Quantifiers

Quantifiers are measures used to characterize specific properties of a system. 
Currently, RMA provides five quantifiers that can be computed or estimated from a microstate distribution. 
Three of them are estimations of classical Recurrence Quantification Analysis (RQA) measures: **recurrence rate**, **determinism**, and **laminarity**. 
One corresponds to an information-theoretic entropy. 
And, the final quantifier is the **disorder** measure, which is defined using the microstate distribution as the basis for applying the method.

In this section, we describe each of these quantifiers and explain how to compute them using the [`RecurrenceMicrostatesAnalysis`](@ref) package.

All quantifiers implemented in the package inherit from [`QuantificationMeasure`](@ref), and their computation is performed using the [`measure`](@ref) function.

```@docs
QuantificationMeasure
measure
```

##  Recurrence microstates entropy

The recurrence microstates entropy (RME) was introduced in 2018 and marks the beginning of the RMA framework [Corso2018Entropy](@cite). It is defined as the Shannon entropy of the RMA distribution:

```math
RME = -\sum_{i = 1}^{2^\sigma} p_i^{(n)} \ln p_i^{(n)},
```

where $n$ is the microstate length, $\sigma$ is the number of recurrence elements constrained within microstate (e.g., $\sigma = n^2$ for square microstates), and $p_i^{(n)}$ denotes the probability of the microstate with decimal representation $i$ in the distribution.

In [`RecurrenceMicrostatesAnalysis`](@ref), the RME is implemented by the [`RecurrenceEntropy`](@ref) struct.

```@docs
RecurrenceEntropy
```

Since the output of the [`distribution`](@ref) function is a [`Probabilities`](@ref) object, the package also supports other information or complexity measures, provided by the package [ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/).

##  Recurrence rate
The recurrence rate quantifies the proportion of black points in a recurrence plot, serving as a measure of the relative density of recurrence points in the sparse recurrence matrix [Webber2015Recurrence](@cite). In standard RQA, it is defined as
```math
RR = \frac{1}{K^2} \sum_{i,j=1}^K R_{i,j}.
```

When estimated using RMA, the recurrence rate is defined as the expected value of the recurrence rate over the microstate distribution:
```math
RR = \sum_{i = 1}^{2^\sigma} p_i^{(n)} RR_i^{(n)},
```
where $RR_i^{(n)}$ denotes the recurrence rate of the $i$-th microstate. For square microstates, this quantity is given by
```math
RR_i^{(n)} = \frac{1}{\sigma} \sum_{m,n=1}^n M_{m,n}^{i, (n)},
```
where $\mathbf M_i^{(n)}$ represents the structure of the microstate identified by the decimal index $i$, with size $n$, for a given motif shape (in this case, square).

In [`RecurrenceMicrostatesAnalysis`](@ref), this quantifier is implemented by the [`RecurrenceRate`](@ref) struct.

```@docs
RecurrenceRate
```

##  Determinism
In standard RQA, determinism (DET) is defined as the fraction of recurrence points that form diagonal line structures in the recurrence plot [Webber2015Recurrence](@cite):
```math
DET = \frac{\sum_{l=d_{min}}^K l~H_D(l)}{\sum_{i,j=1}^N R_{i,j}},
```
where $H_D(l)$ denotes the histogram of diagonal line lengths in the RP, given by
```math
H_D(l)=\sum_{i,j]1}^N(1-R_{i-1,j-1})(1-R_{i+l,j+l})\prod_{k=0}^{l-1}R_{i+k,j+k}.
```

The estimation of DET using RMA is based on the work *"Density-Based Recurrence Measures from Microstates"* [daCruz2025RQAMeasures](@cite). In that work, the DET expression is rewritten as
```math
DET = 1 - \frac{1}{N^2~RR}\sum_{l=1}^{l_{min}-1} l~H_D(l),
```
and the diagonal histogram $H_D(l)$ is related to the RMA distribution through correlations between microstate structures:
```math
\frac{H_D(l)}{(N-l-1)^2}=\vec d^{(l)}\cdot\mathcal{R}^{(l+2)}\vec p^{(l+2)}.
```

For the commonly used case $l_{min} = 2$ (currently the only case implemented in the package), this leads to the approximation
```math
DET\approx 1 - \frac{\vec d^{(1)}\cdot\mathcal{R}^{(3)}\vec p^{3}}{RR}.
```

The correlation term $\vec d^{(1)}\cdot\mathcal{R}^{(3)}\vec p^{3}$ can be simplified by explicitly identifying the microstates selected by $\vec d^{(1)}$.
These correspond to microstates of the form
```math
\begin{pmatrix}
\xi & \xi & 0 \\
\xi & 1 & \xi \\
0 & \xi & \xi \\
\end{pmatrix},
```
where $\xi$ denotes an unconstrained entry. There are $64$ such microstates among the $512$ possible square microstates of size $n = 3$. 
Defining the class $C_D$ as the set of microstates with this structure, DET can be estimated as:
```math
DET\approx 1 - \frac{\sum_{i\in C_D} p_i^{(3)}}{RR},
```
where $p_i^{(3)}$ is the probability of the $i$-th microstate in an RMA distribution of square microstates with size $n = 3$.

A futher simplification can be obtained by defining [`Diagonal`](@ref)-shaped microstates [Ferreira2025RMALib](@cite). 
In the structure above, the unconstrained entries $\xi$ may represent either recurrences or non-recurrences, leading to the need for all $64$ combinations.
Diagonal microstates focus directly on the relevant information, namely the diagonal pattern $0~1~0$.
In this case, DET can be approximated as
```math
DET\approx 1 - \frac{p_3^{(3)}}{RR},
```
where $p_3^{(3)}$ is the probability of observing the diagonal motif $0~1~0$.

In [`RecurrenceMicrostatesAnalysis`](@ref), the computation of DET is implemented by the [`Determinism`](@ref) struct.
```@docs
Determinism
```

##  Laminarity

Laminarity (LAM) is another classical RQA quantifier that measures the proportion of recurrence points forming vertical (line) structures in a recurrence plot. 
It is defined as
```math
LAM = \frac{\sum_{l=v_{min}}^K l~H_V(l)}{\sum_{i,j=1}^N R_{i,j}},
```
where
```math
H_V(l)=\sum_{i,j]1}^N(1-R_{i,j-1})(1-R_{i,j+l})\prod_{k=0}^{l-1}R_{i,j+k}.
```

The estimation of LAM using RMA is also based on the work *"Density-Based Recurrence Measures from Microstates"* [daCruz2025RQAMeasures](@cite) and follows the same logical used for determinsm.
In this case, the estimation requires microstates of the form
```math
\begin{pmatrix}
0 & 1 & 0 \\
\xi & \xi & \xi \\
\xi & \xi & \xi \\
\end{pmatrix},
```
which defines the class $C_L$ of microstates used to estimate LAM as
```math
LAM\approx 1 - \frac{\sum_{i\in C_L} p_i^{(3)}}{RR}.
```

As with determinism, this process can be further simplified by defining a `line` motif [Ferreira2025RMALib](@cite), which captures only the relevant information, namely vertical line patterns of the form $0~1~0$ in the recurrence plot.
In this case, LAM can be approximated as
```math
LAM\approx 1 - \frac{p_3^{(3)}}{RR},
```
where $p_3^{(3)}$ denotes the probability of observing the line motif $0~1~0$.

In [`RecurrenceMicrostatesAnalysis`](@ref), the computation of LAM is implemented by the [`Laminarity`](@ref) struct.

```@docs
Laminarity
```

##  Disorder
The disorder quantifier is implemented based on the work *“Quantifying Disorder in Data”* [Flauzino2025Disorder](@cite). 
It is a novel and powerful tool for data analysis, allowing the differentiation between stochastic and deterministic time series, as well as between different types of stochastic dynamics, such as white, pink, and red Gaussian noise.

Disorder is implemented using square recurrence microstates, which can be permuted by rows and columns and transposed (see [Permutations and Transposition](@ref)). 
This procedure generates a set of equivalent microstates given by
```math
\mathcal{M}_a(\mathbf{M}) = \bigcup_{\sigma_i,\sigma_j\in S_N}\{\mathcal{L}_{\sigma_j}\mathcal{T}\mathcal{L}_{\sigma_i}\mathbf{M},\quad\mathcal{T}\mathcal{L}_{\sigma_j}\mathcal{T}\mathcal{L}_{\sigma_i}\mathbf{M}\}.
```
This defines an equivalence class of microstates denoted by $\mathcal{M}_a$.

The probability of observing a given microstate $\mathbf M_i^{(n)}$ in the recurrence plot, denoted by $p_i^{(n)}$, is computed using [`RecurrenceMicrostatesAnalysis`](@ref).
To compute disorder, the probabilities of microstates belonging to the same class must be normalized.
Thus, for $\mathbf M_i^{(n)} \in \mathcal{M}_a$, the normalized probability within the class is defined as
```math
p_i^{(a, n)} = \frac{p_i^{(n)}}{\sum_{\mathbf{M}_j^{(n)} \in \mathcal{M}_a}~p_j^{(n)}}.
```

The information entropy associated with the probability distribution of microstates in the class $\mathcal{M}_a$ is then defined as
```math
\xi_a(\varepsilon) \stackrel{\mathrm{def}}{=} -\sum_{\mathbf{M}_i^{(n)} \in \mathcal{M}_a} p_i^{(a, n)} \ln p_i^{(a, n)}.
```
This entropy is normalized by $\ln m_a$, where $m_a$ is the number of microstates in the class $\mathcal{M}_a$.
Using [`RecurrenceMicrostatesAnalysis`](@ref), the normalized quantity $\xi_a(\varepsilon) / \ln m_a$ can be computed as
```@example disorder
using Distributions, RecurrenceMicrostatesAnalysis
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
dist = distribution(data, 0.27, 4; sampling = Full())

class = 102
measure(Disorder(4), class, dist)
```

The total entropy over all classes for a given threshold $\varepsilon$ is defined as
```math
\xi(\varepsilon) \stackrel{\mathrm{def}}{=} \frac{1}{A} \sum_{a = 1}^A \frac{\xi_a(\varepsilon)}{\ln m_a},
```
where $A$ is the number of contributing classes and defines the maximum possible amplitude.
This normalization factor can also be computed using [`RecurrenceMicrostatesAnalysis`](@ref):
```@example disorder
A = RecurrenceMicrostatesAnalysis.get_disorder_norm_factor(Disorder(4), data)
```

And the total entropy:
```@example disorder
measure(Disorder(4), dist, A)
```

Finally, the the *disorder index via symmetry in recurrence microstates* (DISREM), or simply **disorder**, is defined as
```math
\Xi = \max_{\varepsilon} \xi(\varepsilon).
```

In [`RecurrenceMicrostatesAnalysis`](@ref), this quantifier is implemented by the [`Disorder`](@ref) struct.

```@docs
Disorder
```

### Computing disorder for compatible time series

Consider a scenario in which a long time series is split into multiple windows. [`RecurrenceMicrostatesAnalysis`](@ref) provides a compact interface to compute the disorder for each window.

As an example, consider a time series with 50,000 points consisting of a sine wave with added white noise, alternating every five windows:
```@example disorder
function data_gen(t)
    x = sin.(6*π .* t)

    count = 0
    for i in 1:1000:50_000
        if count < 5
            x[i:(i-1)+1000] .+= rand(Normal(0, 0.25), 1000)
        elseif count ≥ 9
            count = -1
        end

        count += 1
    end

    return x
end
```
```@example disorder
using CairoMakie

t = range(0, 50, 50_000)
data = data_gen(t)

lines(t, data)
```

The disorder can be computed using the following method:
```julia
measure(settings::Disorder{N}, dataset::Vector{StateSpaceSet}, th_min::Float64, th_max::Float64)
```

To apply it, the time series must first be split into a vector of [`StateSpaceSet`](@ref) objects:
```@example disorder
windows = [ data[(i + 1):(i + 1000)] for i in 0:1000:(length(data) - 1000)]
dataset = Vector{StateSpaceSet}(undef, length(windows))
for i ∈ eachindex(windows)
    dataset[i] = StateSpaceSet(windows[i])
end

dataset
```

Next, the threshold range `th_min` and `th_max` must be defined. There are two possible approaches:
1. Use the full range of admissible threshold values by setting `th_min = 0` and `th_max = maximum(pairwise(Euclidean(), data, data))`, and choosing a small step size via the `num_tests` keywork argument (e.g., `num_tests = 1000`). This approach yields the global maximum disorder values but can be computationally expensive.

2. Use a small interval centered around a known threshold value. This is the recommended approach and is adopted here.

To obtain a suitable reference threshold, we select a subset of windows and compute the optimal disorder threshold using the [`optimize`](@ref) function:
```@example disorder
using Statistics

function find_threshold(disorder, data)
    ths = zeros(Float64, 10)
    for i ∈ eachindex(ths)
        idx = rand(1:length(windows))
        ths[i] = optimize(Threshold(), disorder, data[idx])[1]
    end

    μ = mean(ths)
    σ = std(ths)

    return (max(0.0, μ - 1.5 * σ), μ + 1.5 * σ)
end
```
```@example disorder
dis = Disorder(4)
th_min, th_max = find_threshold(dis, dataset)
```

Finally, the disorder can be computed for all windows using the [`measure`](@ref) function:
```@example disorder
results = measure(dis, dataset, th_min, th_max)
```

```@example disorder
scatterlines(results)
```

!!! tip
    Disorder can also be computed using the [`GPU`](@ref) backend:
    ```julia
    measure(settings::Disorder{N}, dataset::Vector{<:AbstractGPUVector{SVector{D, Float32}}}, th_min::Float32, th_max::Float32)
    ```
    The procedure is identical, but each window must first be transferred to the GPU:
    ```julia
    for i ∈ eachindex(windows)
        dataset[i] = StateSpaceSet(Float32.(windows[i])) |> CuVector
    end
    ```