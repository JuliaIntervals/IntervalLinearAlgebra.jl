module IntervalLinearAlgebra

using StaticArrays, Requires, Reexport
using LinearAlgebra: checksquare

import Base: +, -, *, /, \, ==,
            show, convert, promote_rule, zero, one,
            getindex, IndexStyle, setindex!, size
import CommonSolve: solve

@reexport using LinearAlgebra, IntervalArithmetic

const IA = IntervalArithmetic

export
    set_multiplication_mode, get_multiplication_mode,
    LinearKrawczyk, Jacobi, GaussSeidel, GaussianElimination, HansenBliekRohn, NonLinearOettliPrager, LinearOettliPrager,
    NoPrecondition, InverseMidpoint, InverseDiagonalMidpoint,
    solve, enclose, epsilon_inflation,
    comparison_matrix, interval_norm, interval_isapprox, Orthants,
    is_H_matrix, is_strongly_regular, is_strictly_diagonally_dominant, is_Z_matrix, is_M_matrix,
    rref,
    eigenbox, Rohn, Hertz, verify_eigen, bound_perron_frobenius_eigenvalue,
    AffineExpression, @affinevars,
    AffineParametricArray, AffineParametricMatrix, AffineParametricVector,
    Skalna06


include("linear_systems/enclosures.jl")
include("linear_systems/precondition.jl")
include("linear_systems/solve.jl")
include("linear_systems/verify.jl")
include("linear_systems/oettli.jl")
include("multiplication.jl")
include("utils.jl")
include("classify.jl")
include("rref.jl")
include("pils/affine_expressions.jl")
include("pils/affine_parametric_array.jl")
include("pils/pils_solvers.jl")

include("eigenvalues/interval_eigenvalues.jl")
include("eigenvalues/verify_eigs.jl")

using LinearAlgebra

if Sys.ARCH == :x86_64
    using OpenBLASConsistentFPCSR_jll
end

function  __init__()
    @require IntervalConstraintProgramming = "138f1668-1576-5ad7-91b9-7425abbf3153" include("linear_systems/oettli_nonlinear.jl")
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" include("linear_systems/oettli_linear.jl")
    if Sys.ARCH == :x86_64
        @info "Switching to OpenBLAS with ConsistentFPCSR = 1 flag enabled, guarantees
        coherent floating point rounding mode over all threads"
        BLAS.lbt_forward(OpenBLASConsistentFPCSR_jll.libopenblas_path; verbose =  true)
    else
        BLAS.set_num_threads(1)
        @warn "The number of BLAS threads was set to 1 to ensure rounding mode is consistent"
    end
end

set_multiplication_mode(config[:multiplication])

end
