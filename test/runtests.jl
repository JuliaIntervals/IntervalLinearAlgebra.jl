using IntervalLinearAlgebra, StaticArrays, LazySets
using Test

const IA = IntervalArithmetic

include("test_classify.jl")
include("test_multiplication.jl")
include("test_utils.jl")


include("test_eigenvalues/test_interval_eigenvalues.jl")
include("test_eigenvalues/test_verify_eigs.jl")
include("test_solvers/test_enclosures.jl")
include("test_solvers/test_epsilon_inflation.jl")
include("test_solvers/test_precondition.jl")
include("test_solvers/test_oettli_prager.jl")
