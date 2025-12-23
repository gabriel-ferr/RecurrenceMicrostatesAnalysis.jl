cd(@__DIR__)
include("../src/RecurrenceMicrostatesAnalysis.jl")
using ComplexityMeasures
using Documenter
using DocumenterCitations
using .RecurrenceMicrostatesAnalysis
using StateSpaceSets

pages = [
    "Welcome" => "index.md",
    "Tutorial" => [
        "Distributions" => "tutorial/distributions.md",
        "Quantifiers" => "tutorial/quantifiers.md",
        "Operations" => "tutorial/operations.md",
        "Recurrence Functions" => "tutorial/recurrences.md",
        "Shapes and Sampling" => "tutorial/shapes_and_sampling.md",
        "GPU" => "tutorial/gpu.md",
        "Performance Tips" => "tutorial/tips.md",
        "Utils" => "tutorial/utils.md",
    ],
    "Examples" => [
            "Machine Learning" => "examples/ml.md",
        ],
    "Developers" => "dev.md",
    "References" => "refs.md",
]

bib = CitationBibliography(
    joinpath(@__DIR__, "refs.bib");
    style=:authoryear
)

makedocs(
    sitename = "RecurrenceMicrostatesAnalysis.jl",
    format = Documenter.HTML(
        prettyurls = true,
        collapselevel = 3,
    ),
    authors = "Gabriel Vinicius Ferreira and Felipe Eduardo Lopes da Cruz and Gabriel Marghoti and Thiago de Lima Prado and Sergio Roberto Lopes and Norbert Marwan and Jürgen Kurths",
    modules = [RecurrenceMicrostatesAnalysis, StateSpaceSets, ComplexityMeasures],
    pages = pages,
    doctest = false,
    checkdocs = :exported,
    warnonly = [:doctest, :missing_docs],
    plugins = [bib]
)