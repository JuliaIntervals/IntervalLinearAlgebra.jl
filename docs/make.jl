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
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "wip.md",
        "Applications" => "wip.md",
        "Explanations" => "wip.md",
        "API" => [
            "Interval matrices classification" => "api/classify.md",
            "solver interface" => "api/solve.md",
            "Interval linear systems" => "api/algorithms.md",
            "Preconditioners" => "api/precondition.md",
            "Verified real linear systems" => "api/epsilon_inflation.md",
            "Miscellaneous" => "api/misc.md"
        ],
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/lucaferranti/IntervalLinearAlgebra.jl", devbranch = "main"
)
