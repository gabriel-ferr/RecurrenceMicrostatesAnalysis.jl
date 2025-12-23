export Disorder

##########################################################################################
#   Quantification Measure: Disorder
##########################################################################################
"""
    Disorder{N} <: QuantificationMeasure

Defines the *Disorder* quantifier for microstates of size `N`. The structure contains the `labels`, which identify the microstates belonging to each equivalence class \$\\mathcal{M}_a\$.

To initialize the *Disorder* struct, use:
```julia
Disorder(N)
```
Here, \$N\$ must be equal to 2, 3, 4, or 5. Computing disorder for large values of \$N\$ is currently not supported, as it would require a prohibitive amount of computational memory with the current implementation.

The computation of *Disorder* is performed using the [`measure`](@ref) function:
```julia
measure(settings::Disorder{N}, [x]; kwargs...)
```

### Input
- The `QuantificationMeasure`.
- `[x]`: time-series data provided as a [`StateSpaceSet`](@ref).

### Output
Returns a `Float64` corresponding to the disorder value (\$\\Xi\$).

### Keyword arguments
- `th`: threshold value used to maximize disorder. To improve computational performance, this parameter defines a reference value that limits the search range of thresholds. By default, it is set to the threshold that maximizes disorder for a sampling rate of \$5%\$.
- `th_min`: minimum threshold value defining the search range. By default, this is set to `0.85 * th`.
- `th_max`: maximum threshold value defining the search range. By default, this is set to `1.25 * th`.
- `num_tests`: number of threshold values evaluated within the specified range. The default value is `40`.

### Examples
```julia
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(rand(Uniform(0, 1), 1000))
disrem = measure(Disorder(4), data)
```
"""
struct Disorder{N} <: QuantificationMeasure
    labels::Vector{Vector{Int}}
end

Disorder(N::Int = 3) = 1 < N < 6 ? Disorder{N}(compute_labels(N)) : error("It is not possible to compute disorder for `N = $N`.")
##########################################################################################
#   Implementation: measure
##########################################################################################
function measure(settings::Disorder{N}, class::Int, probs::Probabilities) where {N}
    norm_factor = 0.0
    labels = settings.labels[class]
    @inbounds @simd for i in labels
        norm_factor += probs[i]
    end
    
    if (norm_factor == 0.0)
        return 0.0
    end

    s = 0.0
    @inbounds @simd for i in labels
        p = probs[i] / norm_factor
        s += p * log(p + eps())
    end

    s *= -1
    s /= log(length(labels))

    return s
end

function measure(settings::Disorder{N}, probs::Probabilities, norm_param::Int) where {N}
    total_entropy = 0.0
    for c in 2:length(settings.labels) - 1
        total_entropy += measure(settings, c, probs)
    end

    return total_entropy / norm_param
end

function measure(settings::Disorder{N}, x::StateSpaceSet; th::Float64 = optimize(Threshold(), Disorder(N), x)[1], th_min::Float64 = 0.85 * th, th_max::Float64 = 1.25 * th, num_tests::Int = 40) where {N}
    A = get_disorder_norm_factor(settings, x)
    values = zeros(typeof(th), num_tests)
    th_range = range(th_min, th_max, num_tests)

    for i in eachindex(th_range)
        probs = distribution(x, th_range[i], N; sampling = Full())
        values[i] = measure(settings, probs, A)
    end

    return maximum(values)
end

function measure(settings::Disorder{N}, dataset::Vector{<:AbstractGPUVector{SVector{D, Float32}}}, th_min::Float32, th_max::Float32; num_tests::Int = 40, metric::GPUMetric = GPUEuclidean()) where {N, D}
    A = _norm_factor(Val(N), Val(D))
    values = zeros(Float32, num_tests, length(dataset))
    th_range = Float32.(range(th_min, th_max, num_tests))
    backend = get_backend(dataset[1])

    for i ∈ eachindex(th_range)
        core = GPUCore(backend, Rect(Standard(th_range[i]; metric = metric), N), Full())
        for j in eachindex(dataset)
            probs = distribution(core, dataset[j], dataset[j])
            values[i, j] = measure(settings, probs, A)
        end
    end

    results = zeros(Float32, length(dataset))

    for i ∈ eachindex(dataset)
        results[i] = maximum(values[:, i])
    end

    return results
end

function measure(settings::Disorder{N}, dataset::Vector{StateSpaceSet}, th_min::Float64, th_max::Float64; num_tests::Int = 40, metric::Metric = DEFAULT_METRIC) where {N}
    A = get_disorder_norm_factor(settings, dataset[1])
    values = zeros(Float64, num_tests, length(dataset))
    th_range = range(th_min, th_max, num_tests)

    for i ∈ eachindex(th_range)
        core = CPUCore(Rect(Standard(th_range[i]; metric = metric), N), Full())
        for j in eachindex(dataset)
            probs = distribution(core, dataset[j], dataset[j])
            values[i, j] = measure(settings, probs, A)
        end
    end

    results = zeros(Float32, length(dataset))

    for i ∈ eachindex(dataset)
        results[i] = maximum(values[:, i])
    end

    return results
end
##########################################################################################
#   Utils
##########################################################################################

get_disorder_norm_factor(::Disorder{N}, ::StateSpaceSet{D, T, V}) where {N, D, T, V} = _norm_factor(Val(N), Val(D))
_norm_factor(::Val{2}, ::Val{D}) where D = 4
_norm_factor(::Val{3}, ::Val{D}) where D = D > 1 ? 24 : 23
_norm_factor(::Val{4}, ::Val{D}) where D = D > 1 ? 190 : 145
_norm_factor(::Val{5}, ::Val{D}) where D = D > 1 ? error("Not implemented disorder using N = 5 to data with more than 1 dimension.") : 1173

##########################################################################################
#   Compute labels
##########################################################################################
function compute_labels(N::Int)
    S = collect(permutations(1:N))
    shape = Rect(Standard(0.27), N)

    row_permutation = PermuteRows(shape)
    col_permutation = PermuteColumns(shape; S = S)
    transposition = Transpose(shape)

    identified = Set{Int}()
    labels = Vector{Vector{Int}}(undef, 0)

    for i ∈ 1:(2^(N * N))

        #   Verify if `i` was identified or not.
        if i ∈ identified
            continue
        end

        #   Create a new class.
        push!(labels, Vector{Int}(undef, 1))
        push!(identified, i)
        labels[end][1] = i

        #   Compute the permutations of `i`
        for col in eachindex(col_permutation.Q)
            for row in eachindex(S)

                #   1. Permute rows
                microstate = operate(row_permutation, i, S[row])
                #   2. Permute columns
                microstate = operate(col_permutation, microstate, col)

                #   Security verify
                if microstate ∈ identified
                    continue
                end

                #   Ad the new label to the class
                push!(identified, microstate)
                push!(labels[end], microstate)
            end
        end

        #   Compute transposes
        for label in copy(labels[end])
            
            #   Transpose the label
            microstate = operate(transposition, label)

            #   Security verify
            if microstate ∈ identified
                continue
            end

            #   Ad the new label to the class
            push!(identified, microstate)
            push!(labels[end], microstate)
        end
    end

    return labels
end