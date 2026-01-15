using Test
using RecurrenceMicrostatesAnalysis

@test_throws ArgumentError RecurrenceMicrostatesAnalysis.compute_motif(TestShape(), rand(2, 100), rand(2, 100), [1], [1], SVector{1, Int}(1))
@test_throws ArgumentError RecurrenceMicrostatesAnalysis.get_histogram_size(TestShape())
@test_throws ArgumentError RecurrenceMicrostatesAnalysis.get_power_vector(TestCore(), TestShape())
@test_throws ArgumentError RecurrenceMicrostatesAnalysis.get_offsets(TestCore(), TestShape())