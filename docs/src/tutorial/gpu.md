#   GPU

**RecurrenceMicrostatesAnalysis.jl** supports GPU acceleration for computing RMA distributions.
The GPU backend is implemented using
[KernelAbstractions.jl](https://juliagpu.github.io/KernelAbstractions.jl/stable/), which enables
portable GPU execution across different hardware backends.

The GPU pipeline is implemented via a [`GPUCore`](@ref), which defines a single internal kernel
used to compute microstate histograms across supported devices.

!!! compat
    The GPU kernel is **not compatible with spatial data**.

##  Data requirements
The GPU backend supports **only `Float32` data**. Therefore, input datasets must be explicitly
converted before being used:
```@example gpu
using RecurrenceMicrostatesAnalysis, Distributions
data = StateSpaceSet(Float32.(rand(Uniform(0, 1), 1000)))
```

!!! danger
    When using the GPU backend, inputs must be of type `Float32`. **RecurrenceMicrostatesAnalysis.jl**
    is not compatible with `Float64` on GPU.

##  Recurrence expressions and metrics
When defining a [`RecurrenceExpression`](@ref)  for GPU execution, the distance metric must be a
subtype of [`GPUMetric`](@ref). Metrics from [Distances.jl](https://github.com/JuliaStats/Distances.jl) are
**not supported** on GPU.

For example:
```@example gpu
expr = Standard(0.27f0; metric = GPUEuclidean())
```

!!! compat
    The GPU backend is not compatible with metrics from [Distances.jl](https://github.com/JuliaStats/Distances.jl).

##  Moving data to GPU memory
To enable GPU computation, the data must be transferred to GPU memory. For example:

- Using CUDA:
```julia
using CUDA
gpu_data = data |> CuVector 
```

- Using Metal:
```julia
using Metal
gpu_data = data |> MtlVector
```

###  Output handling
Results computed on the GPU are automatically transferred back to the CPU:
- [`histogram`](@ref) returns a [`Counts`](@ref) object.
- [`distribution`](@ref) returns a [`Probabilities`](@ref) object.

No manual data transfer is required for the output.

##  Metrics for GPU
Since the GPU backend does not support [Distances.jl](https://github.com/JuliaStats/Distances.jl),
distance metrics must be implemented explicitly as subtypes of [`GPUMetric`](@ref).

```@docs
GPUMetric
```

### Implemented GPU metrics
```@docs
GPUEuclidean
```
