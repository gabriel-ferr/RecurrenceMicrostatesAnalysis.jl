cd(@__DIR__)
using ComplexityMeasures
using Documenter
using DocumenterCitations
using RecurrenceMicrostatesAnalysis
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
    modules = [RecurrenceMicrostatesAnalysis, StateSpaceSets, ComplexityMeasures],
    pages = pages,
    doctest = false,
    checkdocs = :exported,
    warnonly = [:doctest, :missing_docs],
    plugins = [bib]
)

deploydocs(
    repo = "github.com/gabriel-ferr/RecurrenceMicrostatesAnalysis.jl.git",
    target = "build",
    push_preview = true
)