#   RMA with Machine Learning

##  Classification using `Flux.jl`
[Flux.jl](https://fluxml.ai/Flux.jl/stable/) is a user-friendly machine learning library in Julia that provides a wide range of tools for building and training neural networks.

In this example, we demonstrate how to use **Flux.jl** together with
**RecurrenceMicrostatesAnalysis.jl** to train a multilayer perceptron (MLP) for classifying
a dynamical system, following the approach presented in [Spezzatto2024ML](@cite).

### Importing required packages
```@example flux
using Flux
using DynamicalSystems
using RecurrenceMicrostatesAnalysis
```

### Generating data
We use the Lorenz system as the data source, which can be simulated using **DynamicalSystems.jl**:
```@example flux
function lorenz!(u, p, t)
    σ, ρ, β = p
    x, y, z = u

    return SVector(
        σ * (y - x),
        x * (ρ - z) - y,
        x * y - β * z
    )
end
```
```@example flux
function lorenz_trajectory(σ, ρ, β; u0 = rand(3), t = 400.0, Ttr = 1200.0, Δt_sample = 0.2)
    p = (σ, ρ, β)
    cds =  ContinuousDynamicalSystem(lorenz!, u0, p)
    x, _ = trajectory(cds, t; Ttr = Ttr, Δt = Δt_sample)
    return x
end
```

We fix the parameters $\sigma = 10$ and $\beta = 8/3$, and vary $\rho \in {26.0, 27.0, 28.0, 29.0, 30.0}$.
The goal is to generate time series for each value of $\rho$ and train a classifier to identify which parameter value generated a given trajectory.

### Creating training and test datasets
First, we define the classes and the size of the training and test datasets:
```@example flux
ρ_cls = [26.0, 27.0, 28.0, 29.0, 30.0];
num_samples_to_test = 50;
num_samples_to_train = 200;
```

First, we define the classes and the size of the training and test datasets:
```@example flux
train_timeseries = Vector{StateSpaceSet}(undef, length(ρ_cls) * num_samples_to_train)
test_timeseries = Vector{StateSpaceSet}(undef, length(ρ_cls) * num_samples_to_test)

train_labels = Vector{Float64}(undef, length(ρ_cls) * num_samples_to_train)
test_labels = Vector{Float64}(undef, length(ρ_cls) * num_samples_to_test)
```

The following function generates the data:
```@example flux
function generate_data!(labels, data, classes, num_elements_per_class)
    c = 1
    for i ∈ eachindex(labels)
        labels[i] = classes[c]
        data[i] = lorenz_trajectory(10.0, classes[c], 8/3)

        if (i % num_elements_per_class == 0)
            c += 1
        end
    end
end
```

```@example flux
generate_data!(train_labels, train_timeseries, ρ_cls, num_samples_to_train)
train_timeseries
```

```@example flux
generate_data!(test_labels, test_timeseries, ρ_cls, num_samples_to_test)
test_timeseries
```

### Preparing the input features
For each time series, we compute the RMA distribution and store it as a feature vector.
```@example flux
microstate_n = 3
train_dataset = Matrix{Float64}(undef, 2^(microstate_n * microstate_n) + 2, length(train_labels))
test_dataset = Matrix{Float64}(undef, 2^(microstate_n * microstate_n) + 2, length(test_labels))
```

The following function computes the RMA features:
```@example flux
function get_probs!(dataset, timeseries, n)
    for i ∈ eachindex(timeseries)
        th, s = optimize(Threshold(), RecurrenceEntropy(), timeseries[i], n)
        dist = distribution(timeseries[i], th, n)
        dataset[1, i] = th
        dataset[2, i] = s
        dataset[3:end, i] = dist[1:end]
    end
end
```

```@example flux
get_probs!(train_dataset, train_timeseries, microstate_n)
train_dataset
```
```@example flux
get_probs!(test_dataset, test_timeseries, microstate_n)
test_dataset
```

### Defining the neural network model
```@example flux
model = Chain(
    Dense(2^(microstate_n * microstate_n) + 2 => 512, identity),
    Dense(512 => 256, selu),
    Dense(256 => 64, selu),
    Dense(64 => length(ρ_cls)),
    softmax
)

model = f64(model)
```

### Training the MLP
First, we encode the labels using one-hot vectors:
```@example flux
train_labels = Flux.onehotbatch(train_labels, ρ_cls)
```
```@example flux
test_labels = Flux.onehotbatch(test_labels, ρ_cls) 
```

We then define the data loader and optimizer:
```@example flux
loader = Flux.DataLoader((train_dataset, train_labels), batchsize = 32, shuffle = true)
```
```@example flux
opt = Flux.setup(Flux.Adam(0.005), model)
```

The training loop is:
```@example flux
for epc ∈ 1:50
    for (x, y) ∈ loader
        _, grads = Flux.withgradient(model) do m
            y_hat = m(x)
            Flux.crossentropy(y_hat, y)
        end

        Flux.update!(opt, model, grads[1])
    end
end
```

### Model evaluation
We compute the classification accuracy as follows:
```@example flux
using LinearAlgebra
function get_quatifiers(predict, trusty, classes)
    conf = zeros(Int, length(classes), length(classes))
    sz = size(predict, 2)

    for i in 1:sz
        mx_prd = findmax(predict[:, i])
        mx_trt = findmax(trusty[:, i])

        conf[mx_prd[2], mx_trt[2]] += 1
    end

    return tr(conf) / (sum(conf) + eps())
end
```
```@example flux
accuracy = get_quatifiers(model(test_dataset), test_labels, ρ_cls)
accuracy * 100
```