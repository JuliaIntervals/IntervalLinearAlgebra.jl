# IntervalLinearAlgebra
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](.LICENSE)
[![Build Status](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/workflows/CI/badge.svg)](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/actions)
[![Coverage](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl/branch/main/graph/badge.svg?token=RYREIXL051)](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/dev)

<p align="center">
    <img src="docs/src/assets/logo.png" alt="IntervalMatrices.jl" width="450"/>
</p>

 <p align="center">
 <i>Linear algebra done rigorously</i></p>

This package contains routines to perform numerical linear algebra using interval arithmetic. At the moment, the package functionalities are limited to interval linear systems, but it is constantly evolving.

## Installation

The package is not registered yet, you can install it from the Julia REPL in package mode as follows

```julia
(@v1.6) pkg> add https://github.com/lucaferranti/intervallinearalgebra.jl
```

## Documentation

Documentation is available [here](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/stable), note that it is still work in progress.

For the time being, have a look at the examples in [examples](./examples/) and [benchmarking](./perf/)

Here is a quick demo about bounding the solution of an interval linear system using a couple of different algorithms.

```julia
julia> using IntervalLinearAlgebra, IntervalArithmetic

julia> A = [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]
4×4 Matrix{Interval{Float64}}:
  [4, 6]   [-1, 1]  [-1, 1]    [-1, 1]
 [-1, 1]  [-6, -4]  [-1, 1]    [-1, 1]
 [-1, 1]   [-1, 1]  [9, 11]    [-1, 1]
 [-1, 1]   [-1, 1]  [-1, 1]  [-11, -9]

julia> b = [-2..4, 1..8, -4..10, 2..12]
4-element Vector{Interval{Float64}}:
  [-2, 4]
   [1, 8]
 [-4, 10]
  [2, 12]

julia> jac = Jacobi()
Jacobi linear solver
max_iterations = 20
atol = 0.0


julia> solve(A, b, jac)
4-element Vector{Interval{Float64}}:
 [-2.60002, 3.10002]
 [-3.90002, 1.65002]
 [-1.48335, 2.15001]
 [-2.35001, 0.794453]

julia>  hbr = HansenBliekRohn()
HansenBliekRohn linear solver


julia> solve(A, b, hbr)
4-element Vector{Interval{Float64}}:
 [-2.50001, 3.10001]
 [-3.90001, 1.20001]
 [-1.40001, 2.15001]
 [-2.35001, 0.6]
```

## References

An excellent introduction to interval linear algebra is
J. Horácek, _Interval Linear and Nonlinear Systems_, 2019, available [here](https://kam.mff.cuni.cz/~horacek/source/horacek_phdthesis.pdf)

See complete the list of [references](./references.md) for the concepts and algorithms used in this package (work in progress).

## Related packages

- [IntervalArithmetic.jl](https://github.com/juliaintervals/IntervalArithmetic.jl) -- Interval computations in Julia
- [IntervalMatrices.jl](https://github.com/JuliaReach/IntervalMatrices.jl) -- Matrices with interval coefficients in Julia.

## Acknowledgment

The development of this package started during the Google Summer of Code (GSoC) 2021 program for the Julia organisation. The author wishes to thank his mentors [David Sanders](https://github.com/dpsanders) and [Marcelo Forets](https://github.com/mforets) for the constant guidance and feedback. During the GSoC program, this project was financially supported by Google.
