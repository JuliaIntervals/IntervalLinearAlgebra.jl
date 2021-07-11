# IntervalLinearAlgebra
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Build Status](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/workflows/CI/badge.svg)](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/actions)
[![Coverage](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl/branch/main/graph/badge.svg?token=RYREIXL051)](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/dev)
<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/stable)-->

<p align="center">
    <img src="docs/src/assets/logo.png" alt="IntervalLinearAlgebra.jl" width="450"/>
</p>

 <p align="center">
 <i>Linear algebra done rigorously</i></p>

## Overview

This package contains routines to perform numerical linear algebra using interval arithmetic. At the moment, the package functionalities are limited to interval linear systems, but it is constantly evolving.

## Installation

The package is not registered yet, you can install it from the Julia REPL in package mode as follows

```julia
(@v1.6) pkg> add https://github.com/lucaferranti/intervallinearalgebra.jl
```

## Documentation
<!-- - [**STABLE**](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/stable) -- Documentation of the latest release -->
- [**DEV**](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/dev) -- Documentation of the current version on main (work in progress)

You can also look at the examples in [examples](./examples/) and [benchmarking](./perf/)

## Quickstart

Here is a quick demo about solving an interval linear system.

```julia
using IntervalLinearAlgebra, IntervalConstraintProgramming, Plots

A = [2..4 -1..1;-1..1 2..4]
b = [-2..2, -1..1]

Xenclose = solve(A, b)
Xexact = solve(A, b, NonLinearOettliPrager())

plot(Xexact.inner, ratio=1, label="exact", legend=:top)
plot!(IntervalBox(Xenclose), label="enclosure")
```
<p align="center">
    <img src="docs/src/assets/quickstart.png" alt="IntervalMatrices.jl" width="450"/>
</p>

## References

An excellent introduction to interval linear algebra is
J. Horácek, _Interval Linear and Nonlinear Systems_, 2019, available [here](https://kam.mff.cuni.cz/~horacek/source/horacek_phdthesis.pdf)

See also the complete list of [references](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/dev/references) for the concepts and algorithms used in this package (work in progress).

## Related packages

- [IntervalArithmetic.jl](https://github.com/juliaintervals/IntervalArithmetic.jl) -- Interval computations in Julia
- [IntervalMatrices.jl](https://github.com/JuliaReach/IntervalMatrices.jl) -- Matrices with interval coefficients in Julia.

## Acknowledgment

The development of this package started during the Google Summer of Code (GSoC) 2021 program for the Julia organisation. The author wishes to thank his mentors [David Sanders](https://github.com/dpsanders) and [Marcelo Forets](https://github.com/mforets) for the constant guidance and feedback. During the GSoC program, this project was financially supported by Google.
