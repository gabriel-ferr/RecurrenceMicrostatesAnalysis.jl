using RecurrenceMicrostatesAnalysis

using Test
using Distances
using Distributions
using RecurrenceAnalysis

const TOLERANCE_RQA = 0.05

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

    @test (abs(det_l2 - measure(Determinism(), x)) / det_l2) ≤ TOLERANCE_RQA

    dist_square = distribution(x, 0.27, 3)
    dist_diagonal = distribution(x, Diagonal(Standard(0.27), 3))

    @test (abs(det_l2 - measure(Determinism(), dist_square)) / det_l2) ≤ TOLERANCE_RQA
    @test (abs(det_l2 - measure(Determinism(), dist_diagonal)) / det_l2) ≤ TOLERANCE_RQA
end

@testset "laminarity" begin
    x = StateSpaceSet(rand(1000))
    rp = RecurrenceMatrix(x, 0.27)
    det_l2 = laminarity(rp)

    @test (abs(det_l2 - measure(Laminarity(), x)) / det_l2) ≤ TOLERANCE_RQA

    dist_square = distribution(x, 0.27, 3)
    dist_diagonal = distribution(x, Diagonal(Standard(0.27), 3))

    @test (abs(det_l2 - measure(Laminarity(), dist_square)) / det_l2) ≤ TOLERANCE_RQA
    @test (abs(det_l2 - measure(Laminarity(), dist_diagonal)) / det_l2) ≤ TOLERANCE_RQA
end

@testset "disorder" begin
    x = StateSpaceSet(rand(Uniform(0, 1), 500))
    @test 0 ≤ measure(Disorder(2), x) ≤ 1
    @test 0 ≤ measure(Disorder(3), x) ≤ 1
    @test 0 ≤ measure(Disorder(4), x) ≤ 1
    # @test 0 ≤ measure(Disorder(5), x) ≤ 1
end