using Documenter, VisualGrowthDynamics

makedocs(;
    modules=[VisualGrowthDynamics],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/skleinbo/VisualGrowthDynamics.jl/blob/{commit}{path}#L{line}",
    sitename="VisualGrowthDynamics.jl",
    authors="Stephan Kleinboelting <skleinbo@thp.uni-koeln.de>",
    assets=String[],
)
