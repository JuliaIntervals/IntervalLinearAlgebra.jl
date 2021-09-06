using IntervalLinearAlgebra
using Plots
using Documenter

DocMeta.setdocmeta!(IntervalLinearAlgebra, :DocTestSetup, :(using IntervalLinearAlgebra); recursive=true)

makedocs(;
    modules=[IntervalLinearAlgebra],
    authors="Luca Ferranti",
    repo="https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="IntervalLinearAlgebra.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaintervals.github.io/IntervalLinearAlgebra.jl",
        assets=String[],
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Linear systems" => "tutorials/linear_systems.md",
            "Eigenvalue computations" => "tutorials/eigenvalues.md"
        ],
        "Applications" => "wip.md",
        "Explanations" => [
            "Interval system solution set" => "explanations/solution_set.md",
            "Preconditioning" => "explanations/preconditioning.md"
        ],
        "API" => [
            "Interval matrices classification" => "api/classify.md",
            "Solver interface" => "api/solve.md",
            "Interval linear systems" => "api/algorithms.md",
            "Preconditioners" => "api/precondition.md",
            "Verified real linear systems" => "api/epsilon_inflation.md",
            "Eigenvalues" => "api/eigenvalues.md",
            "Miscellaneous" => "api/misc.md"
        ],
        "References" => "references.md",
        "Contributing" => "CONTRIBUTING.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaIntervals/IntervalLinearAlgebra.jl", devbranch = "main"
)
