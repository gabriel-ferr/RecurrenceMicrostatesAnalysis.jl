

##########################################################################################
#   Transpose the microstates
##########################################################################################
struct Transpose{R, C} end

Transpose(::Rect2{R, C, B, E}) where {R, C, B, E} = Transpose{R, C}()

##########################################################################################
#   Operate a transposition
##########################################################################################
function operate(::Transpose{R, C}, I::Int) where {R, C}
    B = UInt32(I - 1)
    res_b = UInt32(0)

    for r ∈ 0:R - 1
        for c ∈ 0:C - 1
            res_b |= ((B >> ((r * R) + c)) & 1) << ((c * R) + r)
        end
    end

    return Int(res_b) + 1
end