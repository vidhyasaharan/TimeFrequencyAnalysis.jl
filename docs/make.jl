using TimeFrequencyAnalysis
using Documenter

DocMeta.setdocmeta!(TimeFrequencyAnalysis, :DocTestSetup, :(using TimeFrequencyAnalysis); recursive=true)

makedocs(;
    modules=[TimeFrequencyAnalysis],
    authors="UNSW Sydney, Vidhyasaharan Sethu <v.sethu@unsw.edu.au>",
    sitename="TimeFrequencyAnalysis.jl",
    format=Documenter.HTML(;
        canonical="https://unsw-edu-au.github.io/TimeFrequencyAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/unsw-edu-au/TimeFrequencyAnalysis.jl",
    devbranch="main",
)
