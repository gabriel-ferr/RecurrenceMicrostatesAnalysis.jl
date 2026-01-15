using Test
using RecurrenceMicrostatesAnalysis

@test_throws ArgumentError optimize(Threshold(), RecurrenceEntropy())
@test_throws ArgumentError optimize(Threshold(), Disorder())