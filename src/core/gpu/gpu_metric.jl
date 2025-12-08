export GPUMetric

##########################################################################################
#   GPU Metric
##########################################################################################
#   It is necessary, since Distances.jl doesn't work well with GPU...
abstract type GPUMetric <: Metric end

##########################################################################################
#   Implementation: evaluate
##########################################################################################
function gpu_evaluate(metric::GPUMetric, ::SVector{N, Float32}, ::SVector{N, Float32}) where {N}
    error("Invalid '$(typeof(metric))' implementation.")
end