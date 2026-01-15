using Test
using RecurrenceMicrostatesAnalysis

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname = defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "RecurrenceMicrostatesAnalysis.jl" begin
    #   Core
    testfile("core/backend.jl")
    testfile("core/measures.jl")
    testfile("core/operation.jl")
    testfile("core/optimize.jl")
    testfile("core/recurrence.jl")
    testfile("core/sampling.jl")
    testfile("core/shape.jl")

    #   Recurrences, Sampling, and Shapes
    testfile("recurrences/recurrences.jl")
    testfile("sampling/sampling.jl")
    testfile("shapes/shapes.jl")

    #   RQA
    testfile("rqa/det.jl")
    testfile("rqa/disorder.jl")
    testfile("rqa/entropy.jl")
    testfile("rqa/lam.jl")
    testfile("rqa/rr.jl")

    #   Utils
    testfile("utils/metrics.jl")
    testfile("utils/operations.jl")
    testfile("utils/optimize.jl")

    #   Distributions (all functions)
    testfile("distributions.jl")

    #   GPU Utils
    # testfile("core/gpu_metric.jl")
end