module IntervalLinearAlgebra

# Write your package code here.
using LinearAlgebra, IntervalArithmetic, StaticArrays


export
    Krawczyk, Jacobi, GaussSeidel, GaussElimination, HansenBliekRohn,
    solve, enclose, precondition,
    comparison_matrix


include("solvers/solvers.jl")
include("utils.jl")
end
