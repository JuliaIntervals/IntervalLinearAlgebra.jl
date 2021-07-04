using .LazySets

struct OettliPragerLinear <: DirectSolver end

function (opl::OettliPragerLinear)(A, b)
    n = length(b)
    Ac = mid.(A)
    bc = mid.(b)
    Ar = IntervalArithmetic.radius.(A)
    br = IntervalArithmetic.radius.(b)

    orthants = list_orthants(length(b))
    polyhedra = Vector{HPolyhedron{Float64, Vector{Float64}}}(undef, length(orthants))
    @inbounds for (i, d) in enumerate(orthants)
        D = Diagonal(d)
        Ard = -Ar*D

        A1 = [Ac + Ard; -Ac + Ard; -D]
        b1 = [br + bc; br - bc; zeros(n)]

        polyhedra[i] =  HPolyhedron(A1, b1)
    end
    return polyhedra
end

_default_precondition(_, ::OettliPragerLinear) = NoPrecondition()
