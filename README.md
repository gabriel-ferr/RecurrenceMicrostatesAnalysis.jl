# RecurrenceMicrostatesAnalysis.jl

[![Package Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FRecurrenceMicrostatesAnalysis&query=total_requests&label=Downloads)](http://juliapkgstats.com/pkg/RecurrenceMicrostatesAnalysis)
[![Publication](https://img.shields.io/badge/publication-Chaos-blue.svg)](https://doi.org/10.1063/5.0293708)

**RecurrenceMicrostatesAnalysis.jl** is a simple and fast Julia-based package for recurrence microstates analysis.
It implements the computation of Recurrence Microstates Analysis (RMA) distributions, specific quantifiers — such as [disorder](https://doi.org/10.1103/1y98-x33s) — and the
estimation of typical RQA quantifiers, including determinism and laminarity.

RMA is a subfield of Recurrence Analysis and is a powerful tool for analyzing large time series or large datasets
using statistical methods, offering high performance and avoiding memory issues. Although the field is still
relatively new, it has shown promising applications, including in [Machine Learning](https://doi.org/10.1063/5.0203801).

The package was redesigned in version `0.4.0` to be compatible with [DynamicalSystems](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/)
ecosystem. We therefore recommend exploring the other packages in this ecosystem — expecially [ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/)
and [RecurrenceAnalysis.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/recurrenceanalysis/stable/) —
which can be very useful when working with RMA.

To install the package, run:
```julia
import Pkg
Pkg.add(url="https://github.com/gabriel-ferr/RecurrenceMicrostatesAnalysis.jl")
```

The package documentation is available [online](https://dynamicsufpr.github.io/RecurrenceMicrostatesAnalysis.jl/), or you can build it
locally by running `julia docs/make.jl`.
