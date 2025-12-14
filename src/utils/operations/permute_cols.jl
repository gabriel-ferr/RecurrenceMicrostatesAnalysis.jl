export PermuteColumns

##########################################################################################
#   Permute Rows
##########################################################################################
struct PermuteColumns{R, C}
    Q::Vector{Vector{UInt32}}
end

PermuteColumns(
        ::Rect2{R, C, B, E}; 
        S::Vector{Vector{Int}} = collect(permutations(1:(R * C)))
    ) where {R, C, B, E} = PermuteColumns{R, C}(precompute_Q(R, C, S))

##########################################################################################
#   Operate a permutation of columns
##########################################################################################
function operate(op::PermuteColumns{R, C}, I::Int, Qi::Int) where {R, C}
    B = UInt32(I - 1)
    res_b = UInt32(0)
    mask = UInt32(2^C) - one(UInt32)
    S = op.Q[Qi]

    for n ∈ 0:R - 1
        res_b |= (S[Int((B >> (n * C)) & mask) + 1]) << (n * C)
    end

    return Int(res_b) + 1
end

##########################################################################################
#   Precompute Q
##########################################################################################
function precompute_Q(rows::Int, columns::Int, S::Vector{Vector{Int}})
    Q = Vector{Vector{UInt32}}(undef, length(S))
    B = Vector{UInt32}(undef, rows * columns)

    for i ∈ eachindex(B)
        B[i] = UInt32(i - 1)
    end

    for i ∈ eachindex(S)
        Q[i] = Vector{UInt32}(undef, length(B))
        for n ∈ eachindex(B)
            b = zero(UInt32)
            for c ∈ 0:columns - 1
                b |= ((B[n] >> c) & 1) << (S[i][c + 1] - 1)
            end
            Q[i][n] = b
        end
    end

    return Q
end