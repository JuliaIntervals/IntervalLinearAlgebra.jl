# Linear systems

```@contents
Pages = ["linear_systems.md"]
```

This tutorial will show you how to solve linear systems rigorously using `IntervalLinearAlgebra.jl`.
## Solve interval linear systems

An interval linear system ``\mathbf{Ax}=\mathbf{b}`` is a linear system where ``\mathbf{A}`` and ``\mathbf{b}`` contain intervals. In general, the solution set ``\mathbf{x}`` can have a complex non-convex shape and can thus be hard to characterize exactly (see [this article](../explanations/solution_set.md) for more details). Hence we are interested in finding an interval box containing ``\mathbf{x}``. In `IntervalLinearAlgebra.jl`, this is achieved through the `solve` function, which gives a handy interface to choose the algorithm and preconditioning mechanism. The syntax to call solve is 

```julia
solve(A, b, method, precondition)
```

- ``A`` is an interval matrix
- ``b`` is an interval vector
- `method` is an optional parameter to choose the algorithm used to solve the interval linear system, see below for more details
- `precondition` is an optional parameter to choose the preconditioning for the problem. More details about preconditoining can be found [here](../explanations/preconditioning.md)

### Methods
The supported methods are

- Direct solvers
  - [`GaussianElimination`](@ref)
  - [`HansenBliekRohn`](@ref)
  - [`LinearOettliPrager`](@ref) (requires importing LazySets.jl)

- Iterative solvers
  - [`LinearKrawczyk`](@ref)
  - [`Jacobi`](@ref)
  - [`GaussSeidel`](@ref)
  - [`NonLinearOettliPrager`](@ref) (requires importing IntervalConstraintProgramming.jl)

`LinearOettliPrager` and `NonLinearOettliPrager` are "special" in the sense that they try to exactly characterize the solution set using Oettli-Präger and are not considered in this tutorial. More information about them can be found [here](../explanations/solution_set.md). The other solvers return a vector of intervals, representing an interval enclosure of the solution set. If the method is not specified, Gaussian elimination is used by default.

### Preconditioning

The supported preconditioning mechanisms are

- [`NoPrecondition`](@ref)
- [`InverseMidpoint`](@ref)
- [`InverseDiagonalMidpoint`](@ref)
  
If preconditioning is not specified, then an heuristic strategy based on the type of matrix and solver is used to choose the preconditioning. The strategy is discussed at the end of the [preconditioning tutorial](../explanations/preconditioning.md).

### Examples

We now demonstrate a few examples using the solve function, these examples are taken from [[HOR19]](@ref).

```@example ils
using IntervalLinearAlgebra

A = [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]
```

```@example ils
b = [-2..4, 1..8, -4..10, 2..12]
```

```@example ils
solve(A, b, HansenBliekRohn())
```

```@example ils
solve(A, b, GaussianElimination())
```

```@example ils
solve(A, b, GaussSeidel())
```

For iterative methods, an additional optional parameter `X0` representing an initial guess for the solution's enclosure can be given. If not given, a rough initial enclosure is computed using the [`enclose`](@ref) function.

```@example ils
X0 = fill(-5..5, 4)
solve(A, b, GaussSeidel(), InverseMidpoint(), X0)
```

## Verify real linear systems

`IntervalLinearAlgebra.jl` also offers functionalities to solve real linear systems rigorously. It is of course possible to just convert the real system to an interval system and use the methods described above. In this situation, however, the system will have the property where the diameters of the intervals will be very small (zero or a few floating point units). To solve these kind of systems, it can be more efficient to use the *epsilon inflation* method [[RUM10]](@ref), especially for bigger matrices. Here is an example

```@example ils
A = [1.0 2;3 4]
```

```@example ils
b = [3, 7]
```

the real linear system ``Ax=b`` can now be solved *rigorously* using the [`epsilon_inflation`](@ref) function.

```@example ils
x, cert = epsilon_inflation(A, b)
@show cert
x
```

This function returns two values: an interval vector `x` and a boolean certificate `cert`. If `cert==true` then `x` is guaranteed to be an enclosure of the real linear system `Ax=b`. If `cert == false` then the algorithm could not verify that the enclosure is rigorous, i.e. it may or may not contain the true solution.

In the following example the epsilon inflation method returns a non-rigorous bound

```@example ils
A1 = [1..1+1e-16 2;3 4]
x1, cert = epsilon_inflation(A1, b)
@show cert
x1
```

Since the matrix `A1` is non-regular (it contains the matrix ``\begin{bmatrix}1&2\\3&4\end{bmatrix}`` which is singluar), the solution set is unbounded, hence the algorithm could not prove (rightly) that `x1` is an enclosure of the true solution. 