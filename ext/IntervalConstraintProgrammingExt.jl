module IntervalConstraintProgrammingExt

if !isdefined(Base, :get_extension)
    using ..IntervalLinearAlgebra
    using ..IntervalConstraintProgramming
else
    using IntervalLinearAlgebra
    using IntervalConstraintProgramming
end
using IntervalBoxes: IntervalBox

# Access Symbolics through IntervalConstraintProgramming's dependency
const Symbolics = IntervalConstraintProgramming.Symbolics
const Num = Symbolics.Num

"""
Returns the separator for the Oettli-Prager constraint for the `i`-th row:
`|a_c ⋅x - b_c| - a_r ⋅|x| - b_r <= 0`.

`a_c` and `a_r` are vectors containing midpoints and radii of the intervals in `a`.
Similarly `b_c` and `b_r`.
"""
function oettli_eq(a, b, x)
    ac = mid.(a)
    ar = IntervalArithmetic.radius.(a)
    bc = mid(b)
    br = IntervalArithmetic.radius(b)
    lhs = abs(sum(ac[i] * x[i] for i in eachindex(x)) - bc)
    rhs = sum(ar[i] * abs(x[i]) for i in eachindex(x)) + br
    return constraint(lhs - rhs <= 0)
end

function (op::IntervalLinearAlgebra.NonLinearOettliPrager)(A, b, X::IntervalBox)
    n = length(b)
    vars = [Num(Symbolics.Sym{Real}(Symbol(:x, i))) for i in 1:n]
    separators = [oettli_eq(A[i,:], b[i], vars) for i in 1:n]
    S = reduce(∩, separators)
    return Base.invokelatest(pave, S, X, op.tol)
end

(op::IntervalLinearAlgebra.NonLinearOettliPrager)(A, b, X=enclose(A, b)) = op(A, b, IntervalBox(X...))

IntervalLinearAlgebra._default_precondition(_, ::NonLinearOettliPrager) = NoPrecondition()

end
