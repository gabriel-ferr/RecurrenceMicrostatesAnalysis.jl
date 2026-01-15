using Test
using RecurrenceMicrostatesAnalysis

@test_throws ArgumentError measure(RecurrenceEntropy())
@test_throws ArgumentError measure(RecurrenceRate())
@test_throws ArgumentError measure(Determinism())
@test_throws ArgumentError measure(Laminarity())
@test_throws ArgumentError measure(Disorder())