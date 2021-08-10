# IntervalLinearAlgebra.jl
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/blob/main/LICENSE)
[![Build Status](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/workflows/CI/badge.svg)](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/actions)
[![Coverage](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl/branch/main/graph/badge.svg?token=RYREIXL051)](https://codecov.io/gh/lucaferranti/IntervalLinearAlgebra.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucaferranti.github.io/IntervalLinearAlgebra.jl/dev)

```@raw html
<p align="center">
    <img src="assets/logo.png" alt="IntervalLinearAlgebra.jl" width="450"/>
</p>

 <p align="center">
 <i>Linear algebra done rigorously</i></p>
```

## Overview

Documentation for [IntervalLinearAlgebra.jl](https://github.com/lucaferranti/IntervalLinearAlgebra.jl), a package contains functionalities to solve numerical linear algebra tasks using interval arithmetic.

## Features

!!! note 
    The package is still under active development and everything can change overnight.

- Different algorithms to enclose the solution of an interval linear system
- classify interval matrices
- rigorous solution of real linear systems
- exact characterization of the solution of interval linear systems using Oettli-Präger

## Installation

The package is not registered yet, it can be installed as

```julia
(@v1.6) pkg> add https://github.com/lucaferranti/intervallinearalgebra.jl
```

## Quickstart

```julia
using IntervalLinearAlgebra, LazySets, Plots

A = [2..4 -1..1;-1..1 2..4]
b = [-2..2, -1..1]

Xenclose = solve(A, b)
polytopes = solve(A, b, LinearOettliPrager())

plot(UnionSetArray(polytopes), ratio=1, label="solution set", legend=:top)
plot!(IntervalBox(Xenclose), label="enclosure")
```

![quickstart-example](assets/quickstart.png)
