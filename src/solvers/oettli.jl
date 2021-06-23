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
Converts the interval linear system Ax = b to a system of inequalities using
Oettli-Präger theorem (see section 5.2 of *Jaroslav Horácek, Interval linear and nonlinear systems, 2019*)
and solves the linear system using IntervalConstraintProgramming.jl.

PARAMETERS:
A    : N×N interval matrix

b    : interval vector of length N

X    : Initial bounding box for the solution of Ax = b

vars : list of variables (can be list of symbols, modeling toolkit variables, dynamic polynomaials variables)\n

tol : tolerance of the paving (see documentation of the pave function for further details)
"""
function oettli(A, b, X, vars; tol=0.01)
    separators = [oettli_eq(A[i,:], b[i], vars) for i in 1:length(b)]
    S = reduce(∩, separators)
    return Base.invokelatest(pave, S, X, tol)
end
