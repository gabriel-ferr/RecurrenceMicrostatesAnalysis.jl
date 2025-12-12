using RecurrenceMicrostatesAnalysis

using Distributions
using Metal
using Test

const TOLERANCE = 1e-5

@testset "Metal GPU distributions test" begin
    x = StateSpaceSet(Float32.(rand(1000))) |> MtlVector
    y = StateSpaceSet(Float32.(rand(2000))) |> MtlVector

    @test_nothing distribution(x, y, Rect(Standard(0.27f0; metric = GPUEuclidean()), 2))
    @test_nothing distribution(x, y, Standard(0.27f0; metric = GPUEuclidean()), 2)
    @test_nothing distribution(x, y, 0.27f0, 2)
    @test_nothing distribution(x, Rect(Standard(0.27f0; metric = GPUEuclidean()), 2))
    @test_nothing distribution(x, Standard(0.27f0; metric = GPUEuclidean()), 2)
    @test_nothing distribution(x, 0.27f0, 2)

    dist_1 = distribution(x, 0.27f0, 2)
    dist_2 = distribution(x, 0.27f0, 2)

    @test dist_1 isa Probabilities
    @test dist_2 isa Probabilities

    @test begin
        
        outcomes_1 = outcomes(dist_1)
        outcomes_2 = outcomes(dist_2)

        @test length(outcomes_1) == length(outcomes_2)

        for i ∈ eachindex(outcomes_1)
            if abs(dist_1[i] - dist_2[i]) ≥ TOLERANCE
                return false
            end
        end

        true
    end

    @test_nothing distribution(x, 0.27f0, 2; sampling = Full())
end

@testset "Metal GPU disorder test" begin
    x = rand(Uniform(0, 1), (10000, 3))
    x = Float32.(x)
    v = Vector{MtlVector{SVector{3, Float32}}}(undef, 10)
    for i ∈ 1:10
        v[i] = StateSpaceSet(x[1 + (i - 1) * 1000 : i * 1000, :]) |> MtlVector
    end

    @test begin
        results = measure(Disorder(2), v, 0.25f0, 0.28f0)
        for Ξ ∈ results
            if (abs(1.0 - Ξ) ≥ 0.05)
                return false
            end
        end

        return true
    end 
end