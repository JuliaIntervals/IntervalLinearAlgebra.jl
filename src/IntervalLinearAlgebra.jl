module IntervalLinearAlgebra

using StaticArrays, Requires, Reexport

import CommonSolve: solve

function  __init__()
    @require IntervalConstraintProgramming = "138f1668-1576-5ad7-91b9-7425abbf3153" include("solvers/oettli.jl")
end

@reexport using LinearAlgebra, IntervalArithmetic

export
    Krawczyk, Jacobi, GaussSeidel, GaussElimination, HansenBliekRohn,
    Precondition, NoPrecondition, InverseMidpoint, InverseDiagonalMidpoint,
    solve, enclose, precondition, oettli,
    comparison_matrix, interval_norm, interval_isapprox,
    is_H_matrix, is_strongly_regular, is_strictly_diagonally_dominant, is_Z_matrix, is_M_matrix,
    rref


include("solvers/hull.jl")
include("solvers/precondition.jl")
include("solvers/solve.jl")

include("utils.jl")
include("classify.jl")
include("rref.jl")
end
