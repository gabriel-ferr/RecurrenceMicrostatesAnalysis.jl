export GPUMetric

##########################################################################################
#   GPU Metric
##########################################################################################
#   It is necessary, since Distances.jl doesn't work well with GPU...
"""
    GPUMetric <: Metric

Abstract supertype for metrics compatible with the GPU backend.

Metrics subtyping `GPUMetric` must implement the internal evaluation function
`gpu_evaluate`, which is used during GPU-based computations.

#   Implementations
- [`GPUEuclidean`](@ref)
"""
abstract type GPUMetric <: Metric end

##########################################################################################
#   Implementation: evaluate
##########################################################################################
function gpu_evaluate(metric::GPUMetric, ::SVector{N, Float32}, ::SVector{N, Float32}) where {N}
    error("Invalid '$(typeof(metric))' implementation.")
end

##########################################################################################