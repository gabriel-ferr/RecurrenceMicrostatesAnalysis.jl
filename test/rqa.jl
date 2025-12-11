using RecurrenceMicrostatesAnalysis

using Test
using Distances
using Distributions
using RecurrenceAnalysis

const TOLERANCE = 0.05

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

@testset "determinism" begin
    x = StateSpaceSet(rand(1000))
    rp = RecurrenceMatrix(x, 0.27)
    det_l2 = determinism(rp)

    @test (abs(det_l2 - measure(Determinism(), x)) / det_l2) ≤ TOLERANCE
end

@testset "laminarity" begin
    x = StateSpaceSet(rand(1000))
    rp = RecurrenceMatrix(x, 0.27)
    det_l2 = laminarity(rp)

    @test (abs(det_l2 - measure(Laminarity(), x)) / det_l2) ≤ TOLERANCE
end

@testset "disorder" begin
    x = StateSpaceSet(rand(Uniform(0, 1), 4000))
    @test 0 ≤ measure(Disorder(2), x) ≤ 1
    @test 0 ≤ measure(Disorder(3), x) ≤ 1
    @test 0 ≤ measure(Disorder(4), x) ≤ 1
end