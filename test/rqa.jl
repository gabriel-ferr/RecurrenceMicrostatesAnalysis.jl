using RecurrenceMicrostatesAnalysis

using Test
using Distances

@testset "recurrence rate" begin
    x = StateSpaceSet(rand(1000))
    dist = distribution(x, 0.27, 2)

    @test 0 ≤ measure(RecurrenceRate(), dist) ≤ 1
    @test 0 ≤ measure(RecurrenceRate(), x) ≤ 1
end

@testset "entropy" begin
    x = StateSpaceSet(rand(1000))
    dist = distribution(x, 0.27, 2)

    @test measure(RecurrenceEntropy(), dist) ≥ 0
    @test measure(RecurrenceEntropy(), x) ≥ 0
end