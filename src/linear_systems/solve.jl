"""
    solve(A::AbstractMatrix{T},
          b::AbstractVector{T},
          solver::AbstractIterativeSolver,
          [precondition]::AbstractPrecondition=_default_precondition(A, solver),
          [X]::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

Solves the square interval system ``Ax=b`` using the given algorithm, preconditioner
and initial enclosure

### Input

- `A` -- square interval matrix
- `b` -- interval vector
- `solver` -- algorithm used to solve the linear system
- `precondition` -- preconditioner used. If not given, it is automatically computed based on
                    the matrix `A` and the solver.
- `X` -- initial enclosure.
         if not given, it is automatically computed using [`enclose`](@ref)

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2.0, 4.0]_com  [-1.0, 1.0]_com
 [-1.0, 1.0]_com   [2.0, 4.0]_com

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2.0, 2.0]_com
 [-1.0, 1.0]_com

julia> solve(A, b, GaussSeidel(), NoPrecondition(), [-10..10, -10..10])
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]_trv
 [-1.33334, 1.33334]_trv

julia> solve(A, b, GaussSeidel())
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]_trv_NG
 [-1.33333, 1.33333]_trv_NG
```
"""
function solve(A::AbstractMatrix{T},
               b::AbstractVector{T},
               solver::AbstractIterativeSolver,
               precondition::AbstractPrecondition=_default_precondition(A, solver),
               X::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

    checksquare(A) == length(b) == length(X) || throw(DimensionMismatch())
    A, b = precondition(A, b)
    return solver(A, b, X)

end

"""
    solve(A::AbstractMatrix{T},
          b::AbstractVector{T},
          solver::AbstractDirectSolver,
          [precondition]::AbstractPrecondition=_default_precondition(A, solver)) where
          {T<:Interval}

Solves the square interval system ``Ax=b`` using the given algorithm, preconditioner
and initial enclosure

### Input

- `A` -- square interval matrix
- `b` -- interval vector
- `solver` -- algorithm used to solve the linear system
- `precondition` -- preconditioner used. If not given, it is automatically computed based on
                    the matrix `A` and the solver.

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2.0, 4.0]_com  [-1.0, 1.0]_com
 [-1.0, 1.0]_com   [2.0, 4.0]_com

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2.0, 2.0]_com
 [-1.0, 1.0]_com

julia> solve(A, b, HansenBliekRohn(), InverseMidpoint())
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]_com_NG
 [-1.33333, 1.33333]_com_NG

julia> solve(A, b, HansenBliekRohn())
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]_com
 [-1.33333, 1.33333]_com
```
"""
function solve(A::AbstractMatrix{T},
               b::AbstractVector{T},
               solver::AbstractDirectSolver,
               precondition::AbstractPrecondition=_default_precondition(A, solver)) where
               {T<:Interval}

    checksquare(A) == length(b) || throw(DimensionMismatch())

    A, b = precondition(A, b)
    return solver(A, b)

end


# fallback
"""
    solve(A::AbstractMatrix{T},
          b::AbstractVector{T},
          [solver]::AbstractLinearSolver,
          [precondition]::AbstractPrecondition=_default_precondition(A, solver)) where
          {T<:Interval}

Solves the square interval system ``Ax=b`` using the given algorithm, preconditioner
and initial enclosure

### Input

- `A` -- square interval matrix
- `b` -- interval vector
- `solver` -- algorithm used to solve the linear system. If not given,
              [`GaussianElimination`](@ref) is used.
- `precondition` -- preconditioner used. If not given, it is automatically computed based on
                    the matrix `A` and the solver.

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2.0, 4.0]_com  [-1.0, 1.0]_com
 [-1.0, 1.0]_com   [2.0, 4.0]_com

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2.0, 2.0]_com
 [-1.0, 1.0]_com

julia> solve(A, b)
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]_com
 [-1.33333, 1.33333]_com
```
"""
function solve(A::AbstractMatrix{T},
               b::AbstractVector{T},
               solver::AbstractLinearSolver=_default_solver(),
               precondition::AbstractPrecondition=_default_precondition(A, solver)) where {T<:Interval}

    checksquare(A) == length(b) || throw(DimensionMismatch())

    A, b = precondition(A, b)
    return solver(A, b)

end

## Default settings
_default_solver() = GaussianElimination()

function _default_precondition(A, ::AbstractDirectSolver)
    if is_strictly_diagonally_dominant(A) || is_M_matrix(A)
        return NoPrecondition()
    else
        return InverseMidpoint()
    end
end

# fallback
_default_precondition(_, ::AbstractLinearSolver) = InverseMidpoint()
