using Test
using RecurrenceMicrostatesAnalysis

@test_throws ArgumentError operate(PermuteColumns(Rect(3, 3)))
@test_throws ArgumentError operate(PermuteRows(Rect(3, 3)))
@test_throws ArgumentError operate(Transpose(Rect(3, 3)))