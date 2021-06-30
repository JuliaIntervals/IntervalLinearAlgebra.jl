using IntervalLinearAlgebra
using Documenter

DocMeta.setdocmeta!(IntervalLinearAlgebra, :DocTestSetup, :(using IntervalLinearAlgebra); recursive=true)

makedocs(;
    modules=[IntervalLinearAlgebra],
    authors="Luca Ferranti",
    repo="https://github.com/lucaferranti/IntervalLinearAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="IntervalLinearAlgebra.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucaferranti.github.io/IntervalLinearAlgebra.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "wip.md",
        "Applications" => "wip.md",
        "Explanations" => "wip.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/lucaferranti/IntervalLinearAlgebra.jl", devbranch = "main"
)
