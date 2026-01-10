export GPUEuclidean

##########################################################################################
#   GPUMetric: Euclidean
##########################################################################################
"""
    GPUEuclidean <: GPUMetric

GPU-compatible implementation of the Euclidean distance metric.

```math
d(\\vec{x}, \\vec{y}) = \\sqrt{\\sum_{i = 1}^{m} (x_i - y_i)^2}
```
"""
struct GPUEuclidean <: GPUMetric end

##########################################################################################
#   Metric evaluate
##########################################################################################
@inline function gpu_evaluate(::GPUEuclidean, x, y, i, j, n)
    acc = zero(Float32)

    @inbounds @simd for k in 1:n
        @fastmath d = x[i][k] - y[j][k]
        @fastmath acc += d * d
    end

    return sqrt(acc)
end

##########################################################################################