using IntervalLinearAlgebra, StaticArrays, IntervalConstraintProgramming, LazySets
using Test

const IA = IntervalArithmetic

include("test_classify.jl")
include("test_multiplication.jl")
include("test_utils.jl")

include("test_solvers/test_enclosures.jl")
include("test_solvers/test_epsilon_inflation.jl")
include("test_solvers/test_oettli_prager.jl")
include("test_solvers/test_precondition.jl")
