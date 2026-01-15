using Test
using Distributions
using RecurrenceMicrostatesAnalysis

data = rand(Uniform(0, 1), 1000) |> StateSpaceSet

@test (optimize(Threshold(), RecurrenceEntropy(), data, 2)[1] - 0.27) / 0.27 ≤ 0.1
@test (optimize(Threshold(), Disorder(4), data)[1] - 0.27) / 0.27 ≤ 0.1

@test optimize(Threshold(), RecurrenceEntropy(), data, 2) isa Tuple{Float64, Float64}
@test optimize(Threshold(), Disorder(4), data) isa Tuple{Float64, Float64}