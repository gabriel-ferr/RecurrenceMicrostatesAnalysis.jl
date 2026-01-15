using Test
using Distances
using Distributions
using KernelAbstractions
using StaticArrays
using Random
using RecurrenceMicrostatesAnalysis

@testset "time series" begin
    x = rand(Uniform(0, 1), 1000) |> StateSpaceSet
    y = rand(Uniform(0, 1), 1000) |> StateSpaceSet

    @test sqrt(js_divergence(distribution(x, Rect(Standard(0.27), 3); sampling = SRandom(0.1)), distribution(x, Rect(Standard(0.27), 3); sampling = SRandom(0.1)))) / log(2) ≤ 0.1
    @test sqrt(js_divergence(distribution(x, y, Rect(Standard(0.27), 3); sampling = SRandom(0.1)), distribution(x, y, Rect(Standard(0.27), 3); sampling = SRandom(0.1)))) / log(2) ≤ 0.1
end

@testset "spatial data" begin
    x = rand(Uniform(0, 1), (1, 50, 50))
    y = rand(Uniform(0, 1), (1, 50, 50))

    @test sqrt(js_divergence(distribution(x, Rect(Standard(0.27), (2, 1, 2, 1)); sampling = SRandom(0.05)), distribution(x, Rect(Standard(0.27), (2, 1, 2, 1)); sampling = SRandom(0.05)))) / log(2) ≤ 0.1
    @test sqrt(js_divergence(distribution(x, y, Rect(Standard(0.27), (2, 1, 2, 1)); sampling = SRandom(0.05)), distribution(x, y, Rect(Standard(0.27), (2, 1, 2, 1)); sampling = SRandom(0.05)))) / log(2) ≤ 0.1
end

@testset "GPU" begin
    shape = Rect(Standard(0.27f0), 3)
    core = GPUCore(CPU(), shape, Full())

    @test RecurrenceMicrostatesAnalysis.get_power_vector(core, shape) isa SVector{9, Int32}
    @test RecurrenceMicrostatesAnalysis.get_offsets(core, shape) isa SVector{9, SVector{2, Int32}}
end