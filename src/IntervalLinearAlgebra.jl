module IntervalLinearAlgebra

# Write your package code here.
using LinearAlgebra, IntervalArithmetic, StaticArrays, Requires

function  __init__()
    @require IntervalConstraintProgramming = "138f1668-1576-5ad7-91b9-7425abbf3153" include("solvers/oettli.jl")
end

export
    Krawczyk, Jacobi, GaussSeidel, GaussElimination, HansenBliekRohn,
    solve, enclose, precondition, oettli,
    comparison_matrix, interval_norm, interval_isapprox,
    comparison_matrix, interval_norm


include("solvers/solvers.jl")
include("utils.jl")
end
