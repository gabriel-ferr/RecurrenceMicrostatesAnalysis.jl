using RecurrenceMicrostatesAnalysis
using Test

macro test_nothing(expr)
    return quote
        @test begin
            $(esc(expr))
            true
        end
    end
end

function testfile(file, testname=defaultname(file))
    println("→ running test file $(file)")
    @testset "$testname" begin
        include(file)
    end
end

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))

@testset "RecurrenceMicrostatesAnalysis tests" begin
    testfile("distributions.jl")
    testfile("utils.jl")
    testfile("rqa.jl")

    ##
    ##      GPU:
    # include("gpu/metal.jl")
end