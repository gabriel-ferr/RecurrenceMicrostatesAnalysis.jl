using Test
using Distributions
using RecurrenceAnalysis
using RecurrenceMicrostatesAnalysis

@testset "invalid distribution" begin
    dist = distribution(rand(Uniform(0, 1), 100) |> StateSpaceSet, 0.27, 2)
    @test_throws ArgumentError measure(Determinism(), dist)
end

##  We use a tolerance of 10% here.
@testset "estimation error" begin
    x = StateSpaceSet(rand(Uniform(0, 1), 1000))
    rp = RecurrenceMatrix(x, 0.27)
    det_l2 = determinism(rp)

    @test measure(Determinism(), x) isa Real
    @test (abs(det_l2 - measure(Determinism(), x)) / det_l2) ≤ 0.1

    dist_square = distribution(x, 0.27, 3)
    dist_diagonal = distribution(x, Diagonal(Standard(0.27), 3))

    @test (abs(det_l2 - measure(Determinism(), dist_square)) / det_l2) ≤ 0.1
    @test (abs(det_l2 - measure(Determinism(), dist_diagonal)) / det_l2) ≤ 0.1
end