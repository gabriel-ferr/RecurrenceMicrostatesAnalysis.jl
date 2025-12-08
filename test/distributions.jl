using RecurrenceMicrostatesAnalysis

using Test
using ComplexityMeasures

const TOLERANCE = 1e-5

macro test_nothing(expr)
    return quote
        @test begin
            $(esc(expr))
            true
        end
    end
end

@testset "distributions: StateSpaceSet + CPUCore" begin
    
    x = StateSpaceSet(rand(1000))
    y = StateSpaceSet(rand(2000))

    @test_nothing distribution(x, y, Rect(Standard(0.27), 2))
    @test_nothing distribution(x, y, Standard(0.27), 2)
    @test_nothing distribution(x, y, 0.27, 2)
    @test_nothing distribution(x, Rect(Standard(0.27), 2))
    @test_nothing distribution(x, Standard(0.27), 2)
    @test_nothing distribution(x, 0.27, 2)
    
    dist_1 = distribution(x, 0.27, 2)
    dist_2 = distribution(x, 0.27, 2)

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

    @test_nothing distribution(x, 0.27, 2; sampling = Full())
end

@testset "distributions: AbstractArray + CPUCore" begin
    x = rand(1, 100, 100)
    y = rand(1, 100, 100)

    @test_nothing distribution(x, y, Rect(Standard(0.27), (2, 1, 2, 1)))
    @test_nothing distribution(x, Rect(Standard(0.27), (2, 1, 2, 1)))

    dist_1 = distribution(x, Rect(Standard(0.27), (2, 1, 2, 1)))
    dist_2 = distribution(x, Rect(Standard(0.27), (2, 1, 2, 1)))

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
end