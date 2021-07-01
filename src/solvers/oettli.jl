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
    oettli(A, b, [X]=enclose(A, b); tol=0.01)
Converts the interval linear system ``Ax = b`` to a system of real inequalities which is
solved using `IntervalConstraintProgramming.jl`.

### Input
- `A`   -- N×N interval matrix
- `b`   -- interval vector of length N
- `X`   -- (optional) initial enclosure for the solution of ``Ax = b``. If not given,
            it is automatically computed using [`enclose`](@ref enclose)
- `tol` -- tolerance for paving, default 0.01 (see documentation of `pave` for more details)

### Algorithms

The conversion from interval equalities to real inequalities is based on the Oettli-Präger
theorem [[OET64]](@ref). The resulting system of inequalities is solved using the
Forward-Backward contractor [[JAU14]](@ref).
"""
oettli(A, b, X=enclose(A, b); tol=0.01) = oettli(A, b, IntervalBox(X); tol=tol)

function oettli(A, b, X::IntervalBox; tol=0.01)
    vars = ntuple(i -> Symbol(:x, i), length(b))
    separators = [oettli_eq(A[i,:], b[i], vars) for i in 1:length(b)]
    S = reduce(∩, separators)
    return Base.invokelatest(pave, S, X, tol)
end
