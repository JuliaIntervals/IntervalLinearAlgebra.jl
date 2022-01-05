abstract type ParametricIntervalLinearSolver end

struct Skalna06 <: ParametricIntervalLinearSolver end

function (sk::Skalna06)(A::AffineParametricArray,
                        b::AffineParametricArray,
                        p::Vector{<:Interval})
    mp = map(mid, p)
    R = inv(A(mp))
    bp = b(mp)
    x0 = R * bp
    D = (R * A)(p)

    is_H_matrix(D) || throw(ArgumentError("Could not find an enclosure of given parametric system"))
    z = (R * (b - A * x0))(p)
    Δ = comparison_matrix(D) \ mag.(z)

    return x0 + IntervalArithmetic.Interval.(-Δ, Δ)
end

function solve(A::AffineParametricMatrix,
               b::AffineParametricVector,
               p::Vector{<:Interval},
               solver::ParametricIntervalLinearSolver=Skalna06())

    checksquare(A) == length(b) || throw(DimensionMismatch())

    return solver(A, b, p)
end
