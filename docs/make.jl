using Documenter
using IntervalLinearAlgebra
using Literate
using Plots

const litdir = joinpath(@__DIR__, "literate")

for (root, _, files) in walkdir(litdir)
    for file in files
        if endswith(file, ".jl")
            subfolder = splitpath(root)[end]
            input = joinpath(root, file)
            output = joinpath(@__DIR__, "src", subfolder)
            Literate.markdown(input, output; credit=false,  mdstrings=true)
        end
    end
end
DocMeta.setdocmeta!(IntervalLinearAlgebra, :DocTestSetup, :(using IntervalLinearAlgebra); recursive=true)

makedocs(;
    modules=[IntervalLinearAlgebra],
    authors="Luca Ferranti",
    repo="https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="IntervalLinearAlgebra.jl",
    warnonly=true,
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
        "Applications" => ["Interval FEM" => "applications/FEM_example.md"],
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
    repo="github.com/JuliaIntervals/IntervalLinearAlgebra.jl",
    devbranch = "main",
    push_preview=true
)
