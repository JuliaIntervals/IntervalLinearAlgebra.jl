using .LazySets

function (opl::LinearOettliPrager)(A, b)
    n = length(b)
    Ac = mid.(A)
    bc = mid.(b)
    Ar = IntervalArithmetic.radius.(A)
    br = IntervalArithmetic.radius.(b)

    polytopes = Vector{HPolytope}(undef, 2^n)
    orthants = DiagDirections(n)

    @inbounds for (i, d) in enumerate(orthants)
        D = Diagonal(d)
        Ard = -Ar*D

        A1 = [Ac + Ard; -Ac + Ard; -D]
        b1 = [br + bc; br - bc; zeros(n)]

        polytopes[i] =  HPolytope(A1, b1)
    end
    return identity.(filter!(!isempty, polytopes))
end

_default_precondition(_, ::LinearOettliPrager) = NoPrecondition()
