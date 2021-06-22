module IntervalLinearAlgebra

# Write your package code here.
using LinearAlgebra, IntervalArithmetic, StaticArrays

export
    Krawczyk, Jacobi, GaussSeidel, GaussElimination, HansenBliekRohn,
    solve, enclose, precondition,
    comparison_matrix, interval_norm


include("solvers/solvers.jl")
include("utils.jl")
end
