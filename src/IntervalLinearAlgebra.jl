module IntervalLinearAlgebra

# Write your package code here.
using LinearAlgebra, IntervalArithmetic, StaticArrays, Requires

function  __init__()
    @require IntervalConstraintProgramming = "138f1668-1576-5ad7-91b9-7425abbf3153" include("solvers/oettli.jl")
end

export
    Krawczyk, Jacobi, GaussSeidel, GaussElimination, HansenBliekRohn,
    Precondition, NoPrecondition, InverseMidpointPrecondition, InverseDiagonalMidpointPrecondition,
    solve, enclose, precondition, oettli,
    comparison_matrix, interval_norm, interval_isapprox,
    is_H_matrix, is_strongly_regular, is_strictly_diagonally_dominant, is_Z_matrix, is_M_matrix


include("solvers/solvers.jl")
include("solvers/precondition.jl")
include("utils.jl")
include("classify.jl")
end
