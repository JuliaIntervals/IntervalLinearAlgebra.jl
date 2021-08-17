module IntervalLinearAlgebra

using StaticArrays, Requires, Reexport

import Base: *
import CommonSolve: solve
import IntervalArithmetic: mid
using LinearAlgebra: checksquare

function  __init__()
    @require IntervalConstraintProgramming = "138f1668-1576-5ad7-91b9-7425abbf3153" include("linear_systems/oettli_nonlinear.jl")
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" include("linear_systems/oettli_linear.jl")

    set_multiplication_mode(config[:multiplication])
end

@reexport using LinearAlgebra, IntervalArithmetic

export
    set_multiplication_mode, get_multiplication_mode,
    LinearKrawczyk, Jacobi, GaussSeidel, GaussianElimination, HansenBliekRohn, NonLinearOettliPrager, LinearOettliPrager,
    NoPrecondition, InverseMidpoint, InverseDiagonalMidpoint,
    solve, enclose, epsilon_inflation,
    comparison_matrix, interval_norm, interval_isapprox, list_orthants,
    is_H_matrix, is_strongly_regular, is_strictly_diagonally_dominant, is_Z_matrix, is_M_matrix,
    rref


include("linear_systems/enclosures.jl")
include("linear_systems/precondition.jl")
include("linear_systems/solve.jl")
include("linear_systems/verify.jl")
include("linear_systems/oettli.jl")
include("multiplication.jl")
include("utils.jl")
include("classify.jl")
include("rref.jl")
end
