export PermuteRows

##########################################################################################
#   Permute Rows
##########################################################################################
struct PermuteRows{R, C} end

PermuteRows(::Rect2{R, C, B, E}) where {R, C, B, E} = PermuteRows{R, C}()

##########################################################################################
#   Operate a permutation of rows
##########################################################################################
function operate(::PermuteRows{R, C}, I::Int, σ::Vector{Int}) where {R, C}

    B = UInt32(I - 1)
    res_b = UInt32(0)
    mask = UInt32(2^C) - one(UInt32)

    for n ∈ 0:R - 1
        res_b |= ((B >> (n * C)) & mask) << ((σ[n + 1] - 1) * C)
    end

    return Int(res_b) + 1
end