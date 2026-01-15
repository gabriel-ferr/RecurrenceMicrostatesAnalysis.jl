using Test
using Distributions
using RecurrenceMicrostatesAnalysis

@test_throws ArgumentError Disorder(1)
@test_throws ArgumentError Disorder(6)

@testset "norm factor" begin
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(2), Val(1)) == 4
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(2), Val(2)) == 4
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(3), Val(1)) == 23
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(3), Val(2)) == 24
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(4), Val(1)) == 145
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(4), Val(2)) == 190
    @test RecurrenceMicrostatesAnalysis._norm_factor(Val(5), Val(1)) == 1173
    @test_throws ArgumentError RecurrenceMicrostatesAnalysis._norm_factor(Val(5), Val(2))
end

@testset "basic" begin
    x = StateSpaceSet(rand(Uniform(0, 1), 100))
    @test 0 ≤ measure(Disorder(2), x) ≤ 1
    @test 0 ≤ measure(Disorder(3), x) ≤ 1
    @test 0 ≤ measure(Disorder(4), x) ≤ 1
    # @test 0 ≤ measure(Disorder(5), x) ≤ 1

    @test measure(Disorder(2), x) isa Real
    @test measure(Disorder(3), x) isa Real
    @test measure(Disorder(4), x) isa Real
    # @test measure(Disorder(5), x) isa Real
end

@testset "dataset" begin
    x = rand(Uniform(0, 1), 250)
    windows = [x[(i+1):(i+250)] for i in 0:250:(length(x) - 50)]
    dataset = Vector{StateSpaceSet}(undef, length(windows))
    for i ∈ eachindex(windows)
        dataset[i] = StateSpaceSet(windows[i])
    end

    res = measure(Disorder(2), dataset, 0.25, 0.29)
    for r ∈ res
        @test 0 ≤ r ≤ 1
    end

    @test measure(Disorder(2), dataset, 0.25, 0.29) isa Vector{<:Real}
end