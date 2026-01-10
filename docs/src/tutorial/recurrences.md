#   Recurrence Functions

A recurrence function determines whether two states of a dynamical system,
$\vec{x}$ and $\vec{y}$, are recurrent.

It is defined by an expression of the form
```math
R(\vec{x}, \vec{y}) = \Theta(\varepsilon - |\vec{x} - \vec{y}|),
```
where $\Theta(\cdot)$ denotes the Heaviside step function and $\varepsilon$ is the
threshold parameter defining the maximum distance between two states for them to be
considered $\varepsilon$-recurrent.

This definition differs from the expression used to construct the full recurrence
matrix, introduced in the section [A brief review](@ref). While the recurrence matrix
evaluates all pairwise recurrences in a time series, a recurrence function computes
a single recurrence value between two states.

In **RecurrenceMicrostatesAnalysis.jl**, recurrence functions are implemented via 
the [`RecurrenceExpression`](@ref) abstraction. The actual recurrence evaluation is performed
by implementing the [`recurrence`](@ref) function for the corresponding expression type.

```@docs
RecurrenceExpression
recurrence
```
### Implemented recurrence functions
```@docs
Standard
Corridor
```