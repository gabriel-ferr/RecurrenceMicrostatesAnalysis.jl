export Threshold

##########################################################################################
#   Parameter
##########################################################################################
struct Threshold <: Parameter end

##########################################################################################
#   Implementation: optimize
##########################################################################################
#   - Recurrence Entropy
#.........................................................................................
function optimize(
        ::Threshold, 
        qm::RecurrenceEntropy, 
        x, 
        n::Int; 
        rate::Float64 = 0.05,
        sampling::SamplingMode = SRandom(rate),
        th_max_range::Float64 = 0.5,
        th_start::Float64 = 1e-6,
        fraction::Int = 5,
    )

    ε = th_start
    εopt = 0.0

    if length(x) <= 1000
        εopt = maximum(pairwise(Euclidean(), x)) * (th_max_range - ε)
    else
        εopt = ((maximum(x) - minimum(x)))[1] * size(x, 2)
    end

    Δε = (εopt - ε) / fraction
    fmax = 0.0
    for _ ∈ 1:fraction
        for _ ∈ 1:fraction
            probs = distribution(x, ε, n; sampling = sampling)
            f = measure(qm, probs)

            if f > fmax
                fmax = f
                εopt = ε
            end

            ε += Δε
        end

        ε = εopt - Δε
        Δε *= 2 / fraction
    end

    return εopt, fmax
end

#.........................................................................................
#   - Disorder (random sampling)
#.........................................................................................
function optimize(
        ::Threshold, 
        qm::Disorder{N}, 
        x;
        rate::Float64 = 0.05,
        sampling::SamplingMode = SRandom(rate),
        th_max_range::Float64 = 0.5,
        th_start::Float64 = 1e-6,
        fraction::Int = 5,
    ) where {N}

    ε = th_start
    εopt = 0.0

    if length(x) <= 1000
        εopt = maximum(pairwise(Euclidean(), x)) * (th_max_range - ε)
    else
        εopt = ((maximum(x) - minimum(x)))[1] * size(x, 2)
    end

    A = get_disorder_norm_factor(qm, x)
    Δε = (εopt - ε) / fraction
    fmax = 0.0
    for _ ∈ 1:fraction
        for _ ∈ 1:fraction
            probs = distribution(x, ε, N; sampling = sampling)
            f = measure(qm, probs, A)

            if f > fmax
                fmax = f
                εopt = ε
            end

            ε += Δε
        end

        ε = εopt - Δε
        Δε *= 2 / fraction
    end

    return εopt, fmax
end