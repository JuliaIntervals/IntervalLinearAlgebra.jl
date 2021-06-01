# IntervalLinearAlgebra

[![Build Status](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/workflows/CI/badge.svg)](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/actions)
[![Coverage](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl)

A few experiments of solvers for interval linear systems. So far the following solvers are implemented

- Jacobi
- Gauss-Seidel
- Krawczyk
- Hansen-Bliek-Rohn (HBR)

HBR gives generally tighter intervals and is faster, but it can be used only with H-matrices (or strongly regular matrices using preconditioning in general, but I need to check this).

Example use:

```julia
julia> using IntervalLinearAlgebra, IntervalArithmetic, StaticArrays

julia> A = @SMatrix [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]
4×4 SMatrix{4, 4, Interval{Float64}, 16} with indices SOneTo(4)×SOneTo(4):
  [4, 6]   [-1, 1]  [-1, 1]    [-1, 1]
 [-1, 1]  [-6, -4]  [-1, 1]    [-1, 1]
 [-1, 1]   [-1, 1]  [9, 11]    [-1, 1]
 [-1, 1]   [-1, 1]  [-1, 1]  [-11, -9]

julia> b = @SVector [-2..4, 1..8, -4..10, 2..12]
4-element SVector{4, Interval{Float64}} with indices SOneTo(4):
  [-2, 4]
   [1, 8]
 [-4, 10]
  [2, 12]

julia> jac = Jacobi()
Jacobi(20, 0.0)

julia> solve(A, b, jac)
4-element MVector{4, Interval{Float64}} with indices SOneTo(4):
 [-2.60002, 3.10002]
 [-3.90002, 1.65002]
 [-1.48335, 2.15001]
 [-2.35001, 0.794453]

julia> hbr = HansenBliekRohn()
HansenBliekRohn()

julia> solve(A, b, hbr)
4-element SVector{4, Interval{Float64}} with indices SOneTo(4):
 [-2.50001, 3.10001]
 [-3.90001, 1.2]
 [-1.40001, 2.15001]
 [-2.35001, 0.600001]
```

See the file [`perf/benchmarking.jl`](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/blob/main/perf/benchmark.jl) for some benchmarking. Note that everything is experimental, everything might be good or bad, working or not.

## References

Most algorithms implemented in this package are described in:

[1] J. Horácek, Interval Linear and Nonlinear Systems, 2019, available [here](https://kam.mff.cuni.cz/~horacek/source/horacek_phdthesis.pdf)

## Related packages

- [IntervalMatrices.jl](https://github.com/JuliaReach/IntervalMatrices.jl) -- Matrices with interval coefficients in Julia.

