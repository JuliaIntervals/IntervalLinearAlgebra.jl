@testset "oettli-präger linear method" begin

    A = [2..4 -2..1; -1..2 2..4]
    b = [-2..2, -2..2]

    polyhedra = solve(A, b, LinearOettliPrager())

    for pnt in [[-4, -3], [3, -4], [4, 3], [-3, 4]]
        @test sum(pnt ∈ pol for pol in polyhedra) == 1
    end

    for pnt in [[-5, 5], [5, 5], [5, -5], [-5, -5]]
        @test all(pnt ∉ pol for pol in polyhedra)
    end
end

# NonLinearOettliPrager is disabled: IntervalConstraintProgramming v0.15
# is incompatible with the current Symbolics version.
# @testset "oettli-präger nonlinear method" begin
#     A = [2..4 -2..1; -1..2 2..4]
#     b = [-2..2, -2..2]
#     p = solve(A, b, NonLinearOettliPrager())
#     for pnt in [[-4, -3], [3, -4], [4, 3], [-3, 4]]
#         @test any(pnt ∈ x for x in p.boundary)
#     end
# end
