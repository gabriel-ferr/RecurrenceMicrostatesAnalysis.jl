using Test
using Distances
using Distributions
using KernelAbstractions
using RecurrenceMicrostatesAnalysis

##  Full not implemented for spatial data:
@test_throws ArgumentError distribution(rand(2, 100), Rect(Standard(0.27), (2, 2)); sampling = Full())

@testset "num samples" begin
    data = rand(100) |> StateSpaceSet
    space = SamplingSpace(Rect(2, 2), data, data)

    @test RecurrenceMicrostatesAnalysis.get_num_samples(Full(), space) == space.W * space.H
    @test RecurrenceMicrostatesAnalysis.get_num_samples(Full(), space) isa Integer
end

@testset "get sample" begin
    @testset "CPU" begin
        data_1 = rand(50) |> StateSpaceSet
        data_2 = rand(20) |> StateSpaceSet
        space = SamplingSpace(Rect(2, 2), data_1, data_2)
        core = CPUCore(Rect(Standard(0.27), 3), Full())

        @test RecurrenceMicrostatesAnalysis.get_sample(core, Full(), space, nothing, 10) isa Tuple{<: Integer, <: Integer}
        
        samples = RecurrenceMicrostatesAnalysis.get_num_samples(core.sampling, space)
        for m ∈ 1:samples
            i, j = RecurrenceMicrostatesAnalysis.get_sample(core, Full(), space, nothing, m)
            @test 1 ≤ i ≤ space.W
            @test 1 ≤ j ≤ space.H
        end
    end

    @testset "GPU" begin
        data_1 = rand(50) |> StateSpaceSet
        data_2 = rand(20) |> StateSpaceSet
        space = SamplingSpace(Rect(2, 2), data_1, data_2)
        core = GPUCore(CPU(), Rect(Standard(0.27f0), 3), Full())

        @test RecurrenceMicrostatesAnalysis.get_sample(core, Full(), space, nothing, 10) isa Tuple{<: Integer, <: Integer}

        samples = RecurrenceMicrostatesAnalysis.get_num_samples(core.sampling, space)
        for m ∈ 1:samples
            i, j = RecurrenceMicrostatesAnalysis.get_sample(core, Full(), space, nothing, m)
            @test 1 ≤ i ≤ space.W
            @test 1 ≤ j ≤ space.H
        end
    end
end

@testset "values" begin
    x = StateSpaceSet(rand(100))
    y = StateSpaceSet(rand(50))

    @test js_divergence(distribution(x, 0.27, 3; sampling = Full()), distribution(x, 0.27, 3; sampling = Full())) == 0
    @test js_divergence(distribution(x, y, 0.27, 3; sampling = Full()), distribution(x, y, 0.27, 3; sampling = Full())) == 0
end