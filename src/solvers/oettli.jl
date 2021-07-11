using .IntervalConstraintProgramming

"""
returns the unrolled expression for \$|a ⋅x - b|\$
\$a\$ and \$x\$ must be vectors of the same length and \$b\$ is scalar.
The absolue value in the equation is taken elementwise.
"""
function oettli_lhs(a, b, x)
    ex = :( $(a[1])*$(x[1]))
    for i = 2:length(x)
        ex = :( $ex + $(a[i])*$(x[i]))
    end
    return :(abs($ex - $b))
end

"""
returns the unrolled expression for \$a ⋅|x| + b\$
\$a\$ and \$x\$ must be vectors of the same length and \$b\$ is scalar.
The absolue value in the equation is taken elementwise.
"""
function oettli_rhs(a, b, x)
    ex = :( $(a[1])*abs($(x[1])))
    for i = 2:length(x)
        ex = :( $ex + $(a[i])*abs($(x[i])))
    end
    return :($ex + $b)
end

"""
Returns the separator for the constraint `|a_c ⋅x - b_c| - a_r ⋅|x| - b_r <= 0`.

`a` and `x` must be vectors of the same length and `b` is scalar.

The absolue values in the equation are taken elementwise.

`a_c` and `a_r` are vectors containing midpoints and radii of the intervals in `a`. Similar `b_c` and `b_r`.
"""
function oettli_eq(a, b, x)
    ac = mid.(a)
    ar = radius.(a)
    bc = mid(b)
    br = radius(b)
    lhs = oettli_lhs(ac, bc, x)
    rhs = oettli_rhs(ar, br, x)
    ex = :(@constraint $lhs - $rhs <= 0)
    @eval $ex

end


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
- An object of type `OettliPrager` is a function with methods

        (op::OettliPrager)(A::AbstractMatrix{T},
                           b::AbstractVector{T},
                           [X]::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

        (op::OettliPrager)(A::AbstractMatrix{T},
                           b::AbstractVector{T},
                           X::IntervalBox) where {T<:Interval}

    #### Input
    - `A`   -- N×N interval matrix
    - `b`   -- interval vector of length N
    - `X`   -- (optional) initial enclosure for the solution of ``Ax = b``. If not given,
               it is automatically computed using [`enclose`](@ref enclose)

### Examples

```jldoctest
julia> A = [2..4 -2..1;-1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]

julia> op = NonLinearOettliPrager(0.1)
NonLinearOettliPrager linear solver
tol = 0.1

julia> op(A, b)
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

function (op::NonLinearOettliPrager)(A, b, X::IntervalBox)
    vars = ntuple(i -> Symbol(:x, i), length(b))
    separators = [oettli_eq(A[i,:], b[i], vars) for i in 1:length(b)]
    S = reduce(∩, separators)
    return Base.invokelatest(pave, S, X, op.tol)
end

(op::NonLinearOettliPrager)(A, b, X=enclose(A, b)) = op(A, b, IntervalBox(X))

_default_precondition(_, ::NonLinearOettliPrager) = NoPrecondition()
