#   Utils
This section describes utility functionalities provided by
**RecurrenceMicrostatesAnalysis.jl**, including parameter optimization and operations on
recurrence microstates.

##  Optimizing a parameter
When working with recurrence plots (RPs), a well-known challenge is selecting an appropriate
recurrence threshold. Within the RMA framework, this issue is addressed by optimizing an
information-theoretic or complexity-based measure with respect to the threshold parameter.

In practice, it is common to determine the threshold that maximizes either the
[`RecurrenceEntropy`](@ref) or the [`Disorder`](@ref) quantifier
[Prado2023Sampling](@cite).

This optimization procedure is implemented via the [`optimize`](@ref) function, which
computes the optimal value of a given [`Parameter`](@ref).

```@docs
Parameter
optimize
```

### Parameters
```@docs
Threshold
```

##   Operations on microstates

The package also provides a set of operations that can be applied to recurrence microstates
or to their decimal representations.
```@docs
Operation
operate
```

###  Permutations and Transposition
When working with square microstates, it is natural to consider symmetry operations such as
row and column permutations, as well as transposition. These operations are particularly
important for defining equivalence classes used in the computation of the [`Disorder`](@ref) quantifier.

To illustrate these operations, we consider a square microstate of size $3 \times 3$:
```math
\mathbf{M} = \begin{pmatrix}
a & b & c \\
d & e & f \\
g & h & i
\end{pmatrix}.
```

!!! info
    The primary purpose of these operations is to define equivalence classes of microstates
    used in the computation of [`Disorder`](@ref).

#### Permutations of Rows
Let $\sigma \in S_N$ be a permutation, and let $\mathcal{L}_\sigma$ denote the operator that
permutes the rows of a microstate $\mathbf{M}$ according to $\sigma$.

For example, for $\sigma = 132$, the third and second rows are exchanged, while the first row
remains unchanged:
```math
\mathcal{L}_{132}\mathbf{M} = \begin{pmatrix}
a & b & c \\
g & h & i \\
d & e & f
\end{pmatrix}.
```

For a given microstate size $n$, all possible row (or column) permutations can be generated using the [Combinatorics.jl](https://juliamath.github.io/Combinatorics.jl/stable/) package:
```@example permutation
using Combinatorics

n = 3
Sn = collect(permutations(1:n))
```
```@example permutation
σ = Sn[2]
```

The permutation is applied using [`PermuteRows`](@ref)  operation. For example, consider
the microstate with decimal index $I = 237$:
```math
\mathbf{M} = \begin{pmatrix}
0 & 0 & 1 \\
1 & 0 & 1 \\
1 & 1 & 0
\end{pmatrix}.
```
```@example permutation
using RecurrenceMicrostatesAnalysis

shape = Rect(n, n)
row_permutation = PermuteRows(shape)

operate(row_permutation, 237, σ)
```

The result is the microstate with decimal identifier $I = 239$, corresponding to:
```math
\mathcal{L}_{132}\mathbf{M} = \begin{pmatrix}
0 & 0 & 1 \\
1 & 1 & 0 \\
1 & 0 & 1
\end{pmatrix}.
```

```@docs
PermuteRows
```

#### Permutations of Columns
Column permutations follow the same logic as row permutations, but are applied to the
columns of the microstate.

Let $\mathcal{C}_\sigma$ denote the operator that permutes the columns of $\mathbf{M}$
according to $\sigma \in S_N$. For $\sigma = 132$, the transformation is given by:
```math
\mathcal{C}_{132}\mathbf{M} = \begin{pmatrix}
a & c & b \\
d & f & e \\
g & i & h
\end{pmatrix}.
```

Using the same example of $I = 237$, column permutation is performed using the [`PermuteColumns`](@ref) operation:
```@example permutation
col_permutation = PermuteColumns(shape; S = Sn)
```
```@example permutation
operate(col_permutation, 237, 2)
```

The resulting microstate has decimal identifier $I = 347$, corresponding to:
```math
\mathcal{C}_{132}\mathbf{M} = \begin{pmatrix}
0 & 1 & 0 \\
1 & 1 & 0 \\
1 & 0 & 1
\end{pmatrix}.
```

```@docs
PermuteColumns
```

#### Transposition
Transposition exchanges rows and columns of a microstate. Let $\mathcal{T}$ denote the
transposition operator:
```math
\mathcal{T}\mathbf{M} = \begin{pmatrix}
a & d & g \\
b & e & f \\
c & f & i
\end{pmatrix}.
```

Using the same example microstate with identifier $I = 237$, transposition is performed via the [`Transpose`](@ref) operator:
```@example permutation
transposition = Transpose(shape)
operate(transposition, 237)
```

The resulting microstate has decimal identifier $I = 231$ corresponding to:
```math
\mathcal{T}\mathbf{M} = \begin{pmatrix}
0 & 1 & 1 \\
0 & 0 & 1 \\
1 & 1 & 0
\end{pmatrix}.
```

```@docs
Transpose
```