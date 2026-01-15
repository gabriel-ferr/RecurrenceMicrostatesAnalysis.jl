using Test
using RecurrenceMicrostatesAnalysis

@test operate(PermuteColumns(Rect(3, 3)), 237, 2) == 347
@test operate(PermuteRows(Rect(3, 3)), 237, [1, 3, 2]) == 349
@test operate(Transpose(Rect(3, 3)), 237) == 231