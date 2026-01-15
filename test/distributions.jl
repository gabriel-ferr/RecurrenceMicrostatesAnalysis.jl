using Test
using ComplexityMeasures
using RecurrenceMicrostatesAnalysis

@testset "time series" begin
    x = rand(100) |> StateSpaceSet
    y = rand(200) |> StateSpaceSet
    core = CPUCore(Rect(Standard(0.27), 2), SRandom(0.05))

    @test distribution(core, x, y) isa Probabilities
    @test distribution(core, x) isa Probabilities

    @test distribution(x, 0.27, 3; rate = 0.1) isa Probabilities
    @test distribution(x, 0.27, 3; sampling = Full()) isa Probabilities
    @test distribution(x, 0.27, 3; metric = Cityblock()) isa Probabilities
    @test distribution(x, 0.27, 3; rate = 0.1, metric = Cityblock()) isa Probabilities
    @test distribution(x, 0.27, 3; sampling = Full(), metric = Cityblock()) isa Probabilities

    @test distribution(x, Standard(0.27), 3) isa Probabilities
    @test distribution(x, Standard(0.27), 3; rate = 0.1) isa Probabilities
    @test distribution(x, Standard(0.27), 3; sampling = Full()) isa Probabilities
    @test distribution(x, Standard(0.27; metric = Cityblock()), 3) isa Probabilities
    @test distribution(x, Standard(0.27; metric = Cityblock()), 3; rate = 0.1) isa Probabilities
    @test distribution(x, Standard(0.27; metric = Cityblock()), 3; sampling = Full()) isa Probabilities
    @test distribution(x, Corridor(0.05, 0.27), 3) isa Probabilities
    @test distribution(x, Corridor(0.05, 0.27), 3; rate = 0.1) isa Probabilities
    @test distribution(x, Corridor(0.05, 0.27), 3; sampling = Full()) isa Probabilities
    @test distribution(x, Corridor(0.05, 0.27; metric = Cityblock()), 3) isa Probabilities
    @test distribution(x, Corridor(0.05, 0.27; metric = Cityblock()), 3; rate = 0.1) isa Probabilities
    @test distribution(x, Corridor(0.05, 0.27; metric = Cityblock()), 3; sampling = Full()) isa Probabilities

    @test distribution(x, Rect(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, Rect(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Standard(0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Rect(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Rect(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Standard(0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Triangle(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, Triangle(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Triangle(Standard(0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Triangle(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Triangle(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Triangle(Standard(0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Triangle(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, Triangle(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Triangle(Corridor(0.05, 0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Triangle(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Triangle(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Triangle(Corridor(0.05, 0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities

    @test distribution(x, y, 0.27, 3; rate = 0.1) isa Probabilities
    @test distribution(x, y, 0.27, 3; sampling = Full()) isa Probabilities
    @test distribution(x, y, 0.27, 3; metric = Cityblock()) isa Probabilities
    @test distribution(x, y, 0.27, 3; rate = 0.1, metric = Cityblock()) isa Probabilities
    @test distribution(x, y, 0.27, 3; sampling = Full(), metric = Cityblock()) isa Probabilities

    @test distribution(x, y, Standard(0.27), 3) isa Probabilities
    @test distribution(x, y, Standard(0.27), 3; rate = 0.1) isa Probabilities
    @test distribution(x, y, Standard(0.27), 3; sampling = Full()) isa Probabilities
    @test distribution(x, y, Standard(0.27; metric = Cityblock()), 3) isa Probabilities
    @test distribution(x, y, Standard(0.27; metric = Cityblock()), 3; rate = 0.1) isa Probabilities
    @test distribution(x, y, Standard(0.27; metric = Cityblock()), 3; sampling = Full()) isa Probabilities
    @test distribution(x, y, Corridor(0.05, 0.27), 3) isa Probabilities
    @test distribution(x, y, Corridor(0.05, 0.27), 3; rate = 0.1) isa Probabilities
    @test distribution(x, y, Corridor(0.05, 0.27), 3; sampling = Full()) isa Probabilities
    @test distribution(x, y, Corridor(0.05, 0.27; metric = Cityblock()), 3) isa Probabilities
    @test distribution(x, y, Corridor(0.05, 0.27; metric = Cityblock()), 3; rate = 0.1) isa Probabilities
    @test distribution(x, y, Corridor(0.05, 0.27; metric = Cityblock()), 3; sampling = Full()) isa Probabilities

    @test distribution(x, y, Rect(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, y, Rect(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Rect(Standard(0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Rect(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, y, Rect(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Rect(Standard(0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Triangle(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, y, Triangle(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Triangle(Standard(0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Triangle(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, y, Triangle(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Triangle(Standard(0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Diagonal(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, y, Diagonal(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Diagonal(Standard(0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Diagonal(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, y, Diagonal(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Diagonal(Standard(0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Rect(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, y, Rect(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Rect(Corridor(0.05, 0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Rect(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, y, Rect(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Rect(Corridor(0.05, 0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Triangle(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, y, Triangle(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Triangle(Corridor(0.05, 0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Triangle(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, y, Triangle(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Triangle(Corridor(0.05, 0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Diagonal(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, y, Diagonal(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Diagonal(Corridor(0.05, 0.27), 2); sampling = Full()) isa Probabilities
    @test distribution(x, y, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, y, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, y, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2); sampling = Full()) isa Probabilities
end

@testset "spatial data" begin
    x = rand(1, 20, 20)
    y = rand(1, 10, 10)
    core = CPUCore(Rect(Standard(0.27), (2, 1, 1, 2)), SRandom(0.05))

    @test distribution(x, Rect(Standard(0.27), (2, 1, 1, 2))) isa Probabilities
    @test distribution(x, Rect(Standard(0.27), (2, 1, 1, 2)); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Standard(0.27; metric = Cityblock()), (2, 1, 1, 2))) isa Probabilities
    @test distribution(x, Rect(Standard(0.27; metric = Cityblock()), (2, 1, 1, 2)); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27), 2)) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Diagonal(Standard(0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27), (2, 1, 1, 2))) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27), (2, 1, 1, 2)); rate = 0.1) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27; metric = Cityblock()), (2, 1, 1, 2))) isa Probabilities
    @test distribution(x, Rect(Corridor(0.05, 0.27; metric = Cityblock()), (2, 1, 1, 2)); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27), 2)) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27), 2); rate = 0.1) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2)) isa Probabilities
    @test distribution(x, Diagonal(Corridor(0.05, 0.27; metric = Cityblock()), 2); rate = 0.1) isa Probabilities
end