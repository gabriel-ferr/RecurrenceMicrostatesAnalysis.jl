module RecurrenceMicrostatesAnalysis

@doc let 
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end RecurrenceMicrostatesAnalysis

##########################################################################################
#   Packages and constants
##########################################################################################
using ComplexityMeasures
using Distances
using GPUArraysCore
using Reexport
using StaticArrays

@reexport using StateSpaceSets

const DEFAULT_METRIC = Euclidean()

##########################################################################################
#   Core API types and functions
##########################################################################################
include("core/abstract_core.jl")
include("core/recurrence.jl")
include("core/shape.jl")
include("core/sampling.jl")

include("core/cpu_core.jl")

##########################################################################################
#   Recurrence functions, motif shapes, and sampling modes
##########################################################################################
include("recurrences/recurrences.jl")
include("shapes/shapes.jl")
include("sampling/sampling.jl")

##########################################################################################
#   Quantifiers
##########################################################################################

##########################################################################################
#   Utils
##########################################################################################

end