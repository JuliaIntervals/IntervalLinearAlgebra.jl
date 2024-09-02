using IntervalLinearAlgebra, StaticArrays, LazySets, IntervalConstraintProgramming
using Test

const IA = IntervalArithmetic

include("test_classify.jl")
include("test_multiplication.jl")
include("test_utils.jl")

include("test_numerical_test/test_numerical_test.jl")
include("test_eigenvalues/test_interval_eigenvalues.jl")
include("test_eigenvalues/test_verify_eigs.jl")
include("test_eigenvalues/test_miyajima_svd.jl")
include("test_solvers/test_enclosures.jl")
include("test_solvers/test_epsilon_inflation.jl")
include("test_solvers/test_precondition.jl")
include("test_solvers/test_oettli_prager.jl")
include("test_pils/test_linexpr.jl")
include("test_pils/test_affine_parametic_array.jl")
include("test_pils/test_pils_solvers.jl")
