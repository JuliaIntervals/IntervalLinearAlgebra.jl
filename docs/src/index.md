![](assets/logo-text.svg)

![version](https://img.shields.io/github/v/release/juliaintervals/IntervalLinearAlgebra.jl)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/lucaferranti/IntervalLinearAlgebra.jl/blob/main/LICENSE)[![Build Status](https://github.com/juliaintervals/IntervalLinearAlgebra.jl/workflows/CI/badge.svg)](https://github.com/juliaintervals/IntervalLinearAlgebra.jl/actions)[![Coverage](https://codecov.io/gh/juliaintervals/IntervalLinearAlgebra.jl/branch/main/graph/badge.svg?token=mgCzKMPiwK)](https://codecov.io/gh/juliaintervals/IntervalLinearAlgebra.jl)[![bibtex citation](https://img.shields.io/badge/bibtex-citation-green)](#Citation)[![zenodo doi](https://img.shields.io/badge/zenodo-DOI-blue)](https://doi.org/10.5281/zenodo.5363563)

## Overview

This package contains routines to perform numerical linear algebra using interval arithmetic. This can be used both for rigorous computations and uncertainty propagation.

An first overview of the package was given at JuliaCon 2021, the slides are available [here](https://github.com/lucaferranti/ILAjuliacon2021).

```@raw html
<iframe width="560" height="315" src="https://www.youtube.com/embed/fre0TKgLJwg" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

## Features

!!! note 
    The package is still under active development and things evolve quickly (or at least should)

- enclosure of the solution of interval linear systems
- exact characterization of the solution set of interval linear systems using Oettli-Präger
- verified solution of floating point linear systems
- enclosure of eigenvalues of interval matrices
- verified computation of eigenvalues and eigenvectors of floating point matrices

## Installation

Open a Julia session and enter

```julia
using Pkg; Pkg.add("IntervalLinearAlgebra")
```

this will download the package and all the necessary dependencies for you. Next you can import the package with

```julia
using IntervalLinearAlgebra
```

and you are ready to go.

## Quickstart

```julia
using IntervalLinearAlgebra, LazySets, Plots

A = [2..4 -1..1; -1..1 2..4]
b = [-2..2, -1..1]

Xenclose = solve(A, b)
polytopes = solve(A, b, LinearOettliPrager())

plot(UnionSetArray(polytopes), ratio=1, label="solution set", legend=:top)
plot!(IntervalBox(Xenclose), label="enclosure")
```

![quickstart-example](assets/quickstart.png)

## Citation

If you use this package in your work, please cite it as

```
@software{ferranti2021interval,
author = {
            Luca Feranti and
            Marcelo Forets and
            David P. Sanders
         },
title  = {IntervalLinearAlgebra.jl: linear algebra done rigorously},
month  = {9},
year   = {2021},
doi    = {10.5281/zenodo.5363563},
url    = {https://github.com/juliaintervals/IntervalLinearAlgebra.jl}
}
```
