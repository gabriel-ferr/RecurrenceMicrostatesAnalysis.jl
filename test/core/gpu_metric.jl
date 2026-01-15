using Test
using RecurrenceMicrostatesAnalysis

struct TestGPUMetric <: GPUMetric end
data = SVector{1, Float32}(0.0f0)

@test_throws ArgumentError RecurrenceMicrostatesAnalysis.gpu_evaluate(TestGPUMetric(), data, data)