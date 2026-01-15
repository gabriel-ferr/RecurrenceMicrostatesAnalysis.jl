using Test
using RecurrenceMicrostatesAnalysis

@test_throws ArgumentError Standard(-1.0)

@testset "time series" begin
    x = rand(100) |> StateSpaceSet
    y = rand(100) |> StateSpaceSet
    
    standard = Standard(0.27)

    @test RecurrenceMicrostatesAnalysis.recurrence(standard, x, y, 1, 1) isa UInt8

    r = RecurrenceMicrostatesAnalysis.recurrence(standard, x, y, 1, 1)
    @test r == UInt8(0) || r == UInt8(1)
end

@testset "GPU (running on CPU)" begin
    x = rand(100) |> StateSpaceSet
    y = rand(100) |> StateSpaceSet
    
    standard = Standard(0.27; metric = GPUEuclidean())

    @test RecurrenceMicrostatesAnalysis.gpu_recurrence(standard, x, y, 1, 1, 1) isa UInt8

    r = RecurrenceMicrostatesAnalysis.gpu_recurrence(standard, x, y, 1, 1, 1)
    @test r == UInt8(0) || r == UInt8(1)
end

@testset "spatial data" begin
    x = rand(2, 100)
    y = rand(2, 100)
    
    standard = Standard(0.27)

    pos = (1, )
    @test RecurrenceMicrostatesAnalysis.recurrence(standard, x, y, pos, pos) isa UInt8

    r = RecurrenceMicrostatesAnalysis.recurrence(standard, x, y, pos, pos)
    @test r == UInt8(0) || r == UInt8(1)
end