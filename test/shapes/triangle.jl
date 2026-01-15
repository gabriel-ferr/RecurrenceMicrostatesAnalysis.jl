using Test
using Distances
using Distributions
using KernelAbstractions
using StaticArrays
using Random
using RecurrenceMicrostatesAnalysis

##  Triangle not implemented for spatial data:
@test_throws ArgumentError distribution(rand(2, 100), Triangle(Standard(0.27), 3))

@testset "time series" begin
    x = rand(Uniform(0, 1), 1000) |> StateSpaceSet
    y = rand(Uniform(0, 1), 1000) |> StateSpaceSet

    @test sqrt(js_divergence(distribution(x, Triangle(Standard(0.27), 3); sampling = SRandom(0.1)), distribution(x, Triangle(Standard(0.27), 3); sampling = SRandom(0.1)))) / log(2) ≤ 0.1
    @test sqrt(js_divergence(distribution(x, y, Triangle(Standard(0.27), 3); sampling = SRandom(0.1)), distribution(x, y, Triangle(Standard(0.27), 3); sampling = SRandom(0.1)))) / log(2) ≤ 0.1
end

@testset "GPU" begin
    shape = Triangle(Standard(0.27f0), 3)
    core = GPUCore(CPU(), shape, Full())

    @test RecurrenceMicrostatesAnalysis.get_power_vector(core, shape) isa SVector{6, Int32}
    @test RecurrenceMicrostatesAnalysis.get_offsets(core, shape) isa SVector{6, SVector{2, Int32}}
end