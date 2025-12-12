using RecurrenceMicrostatesAnalysis

using Test
using Distributions

const TOLERANCE = 0.15

@testset "optimize threshold" begin
    x = StateSpaceSet(rand(Uniform(0, 1), 2000))

    @test (abs(optimize(Threshold(), RecurrenceEntropy(), x, 3)[1] - 0.27) / 0.27) ≤ TOLERANCE
    @test (abs(optimize(Threshold(), Disorder(), x)[1] - 0.27) / 0.27) ≤ TOLERANCE
end