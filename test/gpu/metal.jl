using RecurrenceMicrostatesAnalysis

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

        for i in eachindex(outcomes_1)
            if abs(dist_1[i] - dist_2[i]) ≥ TOLERANCE
                return false
            end
        end

        true
    end

    @test_nothing distribution(x, 0.27f0, 2; sampling = Full())
end