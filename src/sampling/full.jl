export Full

##########################################################################################
#   Sampling Mode: Random
##########################################################################################
"""
    Full <: SamplingMode

Sampling mode that selects **all** possible microstates within the
[`SamplingSpace`](@ref).

#   Constructor
```julia
s = Full()
```

!!! warning
    The **Full** sampling mode is not supported for spatial data.
"""
struct Full <: SamplingMode end

##########################################################################################
#   Implementation: sampling
##########################################################################################
#   Based on time series: (CPU)
#.........................................................................................
function get_sample(::CPUCore, ::Full, space::SSRect2, _, m)
    i = ((m - 1) % space.W) + 1
    j = ((m - 1) ÷ space.W) + 1

    return i, j
end
#.........................................................................................
#   Based on time series: (GPU)
#.........................................................................................
function get_sample(::GPUCore, ::Full, space::SSRect2, _, m)
    i = Int32((m - 1) % space.W) + 1
    j = Int32((m - 1) ÷ space.W) + 1

    return i, j
end

##########################################################################################
#   Implementations: Utils
##########################################################################################
get_num_samples(::Full, space::SSRect2) = space.W * space.H

##########################################################################################