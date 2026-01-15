using Test
using RecurrenceMicrostatesAnalysis

data = rand(10) |> StateSpaceSet
@test RecurrenceMicrostatesAnalysis.gpu_evaluate(GPUEuclidean(), data, data, 1, 1, 1) == 0