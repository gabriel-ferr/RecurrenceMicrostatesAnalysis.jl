using RecurrenceMicrostatesAnalysis
using Test

function testfile(file, testname=defaultname(file))
    println("→ running test file $(file)")
    @testset "$testname" begin
        include(file)
    end
end

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))

@testset "RecurrenceMicrostatesAnalysis tests" begin
    testfile("distributions.jl")
end