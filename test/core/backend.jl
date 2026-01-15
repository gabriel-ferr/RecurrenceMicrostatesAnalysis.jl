using Test
using RecurrenceMicrostatesAnalysis

struct TestCore <: RMACore end
x = rand(10) |> StateSpaceSet

@test_throws ArgumentError histogram(TestCore(), x, x)
@test_throws ArgumentError distribution(TestCore(), x, x)