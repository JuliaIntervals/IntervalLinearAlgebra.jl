module IntervalLinearAlgebra

# Write your package code here.
using LinearAlgebra, IntervalArithmetic, StaticArrays, IntervalConstraintProgramming

export
    Krawczyk, Jacobi, GaussSeidel, GaussElimination, HansenBliekRohn,
    solve, enclose, precondition, oettli,
    comparison_matrix, interval_norm, interval_isapprox,
    comparison_matrix, interval_norm


include("solvers/solvers.jl")
include("utils.jl")
include("solvers/oettli.jl")
end
