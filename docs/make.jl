using TimeFrequencyAnalysis
using Documenter

DocMeta.setdocmeta!(TimeFrequencyAnalysis, :DocTestSetup, :(using TimeFrequencyAnalysis); recursive=true)

makedocs(;
    modules=[TimeFrequencyAnalysis],
    authors="UNSW Sydney, Vidhyasaharan Sethu <v.sethu@unsw.edu.au>",
    sitename="TimeFrequencyAnalysis.jl",
    format=Documenter.HTML(;
        canonical="https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Data Structures" => "types.md",
        "Framing" => "framing.md",
        "Spectral Analyses" => "spectral.md",
        "Utilities" => "utilities.md",
        "Internals" => "internals.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/vidhyasaharan/TimeFrequencyAnalysis.jl",
    devbranch="main",
)
