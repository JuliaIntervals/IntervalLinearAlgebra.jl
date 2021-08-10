"""
    LinearOettliPrager <: AbstractDirectSolver

Type for the OettliPrager solver of the interval linear system ``Ax=b``. The solver first
converts the system of interval equalities into a system of real inequalities using
Oettli-Präger theorem [[OET64]](@ref) and then finds the feasible set by solving a LP
problem in each orthant using `LazySets.jl`.

### Notes
- You need to import `LazySets.jl` to use this functionality.
- An object of type `LinearOettliPrager` is a function with methods

        (op::LinearOettliPrager)(A::AbstractMatrix{T},
                                 b::AbstractVector{T}) where {T<:Interval}

    #### Input
    - `A`   -- N×N interval matrix
    - `b`   -- interval vector of length N

### Examples

```julia-repl
julia> A = [2..4 -2..1;-1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]

julia> polytopes = solve(A, b, LinearOettliPrager());

julia> typeof(ans)
Vector{HPolytope{Float64, SparseArrays.SparseVector{Float64, Int64}}}
```
"""
struct LinearOettliPrager <: AbstractDirectSolver end


"""
    NonLinearOettliPrager <: AbstractIterativeSolver

Type for the OettliPrager solver of the interval linear system ``Ax=b``. The solver first
converts the system of interval equalities into a system of real inequalities using
Oettli-Präger theorem [[OET64]](@ref) and then finds the feasible set using
the forward-backward contractor method [[JAU14]](@ref) implemented in
`IntervalConstraintProgramming.jl`.

### Fields
- `tol` -- tolerance for the paving, default 0.01.

### Notes
- You need to import `IntervalConstraintProgramming.jl` to use this functionality.
- An object of type `NonLinearOettliPrager` is a function with methods

        (op::NonLinearOettliPrager)(A::AbstractMatrix{T},
                                    b::AbstractVector{T},
                                    [X]::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

        (op::NonLinearOettliPrager)(A::AbstractMatrix{T},
                                    b::AbstractVector{T},
                                    X::IntervalBox) where {T<:Interval}

    #### Input
    - `A`   -- N×N interval matrix
    - `b`   -- interval vector of length N
    - `X`   -- (optional) initial enclosure for the solution of ``Ax = b``. If not given,
               it is automatically computed using [`enclose`](@ref enclose)

### Examples

```julia-repl
julia> A = [2..4 -2..1;-1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]

julia> solve(A, b, NonLinearOettliPrager(0.1))
Paving:
- tolerance ϵ = 0.1
- inner approx. of length 1195
- boundary approx. of length 823
```
"""
struct NonLinearOettliPrager <: AbstractIterativeSolver
    tol::Float64
end
NonLinearOettliPrager() = NonLinearOettliPrager(0.01)
