using IntervalLinearAlgebra, StaticArrays, IntervalConstraintProgramming, LazySets
using Test

@testset "IntervalLinearAlgebra.jl" begin
    @testset "precondition" begin
        A = [2..4 -2..1; -1..2 2..4]
        b = [-2..2, -2..2]

        np = NoPrecondition()
        idp = InverseDiagonalMidpoint()
        imp = InverseMidpoint()

        A1, b1 = np(A, b)
        @test A1 == A && b1 == b

        A2, b2 = idp(A, b)
        Acorrect = [2/3..4/3 -2/3..1/3; -1/3..2/3 2/3..4/3]
        bcorrect = [-2/3..2/3, -2/3..2/3]
        @test all(interval_isapprox.(A2, Acorrect)) && all(interval_isapprox.(b2, bcorrect))

        A3, b3 = imp(A, b)
        Acorrect = [22/37..52/37 -20/37..20/37;-20/37..20/37 22/37..52/37]
        bcorrect = [-28/37..28/37, -28/37..28/37]
        @test all(interval_isapprox.(A3, Acorrect)) && all(interval_isapprox.(b3, bcorrect))
    end

    @testset "Test linear solvers" begin
        As = @SMatrix [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]
        bs = @SVector [-2..4, 1..8, -4..10, 2..12]

        Am = Matrix(As)
        bm = Vector(bs)
        jac = Jacobi()
        gs = GaussSeidel()
        hbr = HansenBliekRohn()
        kra = LinearKrawczyk()

        for (A, b) in zip([As, Am], [bs, bm])
            xgs = solve(A, b, gs)
            xjac = solve(A, b, jac)
            xhbr = solve(A, b, hbr)
            xkra = solve(A, b, kra)

            @test all(interval_isapprox.(xgs, [-2.6..3.1, -3.9..1.65, -1.48..2.15, -2.35..0.79]; atol=0.01))
            @test all(interval_isapprox.(xjac, [-2.6..3.1, -3.9..1.65, -1.48..2.15, -2.35..0.79]; atol=0.01))

            @test all(interval_isapprox.(xhbr, [-2.5..3.1, -3.9..1.2, -1.4..2.15, -2.35..0.6]; atol=0.01))

            @test all(interval_isapprox.(xkra, [-8..8, -8..8, -8..8, -8..8]; atol=0.01))
        end

        ge = GaussianElimination()
        xge = solve(Am, bm, ge)
        @test all(interval_isapprox.(xge, [-2.6..3.1, -3.9..1.5, -1.43..2.15, -2.35..0.6]; atol=0.01))

        xdef = solve(Am, bm)
        @test all(interval_isapprox.(xdef, [-2.6..3.1, -3.9..1.5, -1.43..2.15, -2.35..0.6]; atol=0.01))

        A = [2..4 -2..1; -1..2 2..4]
        b = [-2..2, -2..2]

        x1 = solve(A, b)
        @test all(interval_isapprox.(x1, [-14..14, -14..14]))
        x2 = solve(A, b, HansenBliekRohn())
        @test all(interval_isapprox.(x2, [-14..14, -14..14]))
    end

    @testset "oettli-präger method" begin

        A = [2..4 -2..1; -1..2 2..4]
        b = [-2..2, -2..2]

        p = solve(A, b, NonLinearOettliPrager())

        polyhedra = solve(A, b, LinearOettliPrager())

        for pnt in [[-4, -3], [3, -4], [4, 3], [-3, 4]]
            @test any(pnt ∈ x for x in p.boundary)
            @test sum(pnt ∈ pol for pol in polyhedra) == 1
        end

        for pnt in [[-5, 5], [5, 5], [5, -5], [-5, -5]]
            @test all(pnt ∉ pol for pol in polyhedra)
        end
    end

    @testset "classify matrices" begin
        A = [0..2 1..1;-1.. -1 0..2]
        @test !is_H_matrix(A)
        @test !is_strongly_regular(A)

        B = [-2.. -2 1..1; 5..6 -2.. -2]
        @test is_strongly_regular(B)
        @test !is_Z_matrix(B)
        @test !is_M_matrix(B)

        C = [2..2 1..1; 0..2 2..2]
        @test is_H_matrix(C)
        @test !is_strictly_diagonally_dominant(C)

        D = [2..2 -1..0; -1..0 2..2]
        @test is_strictly_diagonally_dominant(D)
        @test is_Z_matrix(D)
        @test is_M_matrix(D)

        E = [2..4 -2..1;-1..2 2..4]
        @test !is_Z_matrix(E)
        @test !is_M_matrix(E)
        @test !is_H_matrix(E)
        @test is_strongly_regular(E)
    end

    @testset "Gaussian elimination" begin
        A1 = [1..2 1..2;2..2 3..3]
        @test rref(A1) == [2..2 3..3; 0..0 -2..0.5]

        A2 = fill(0..0, 2, 2)
        @test_throws ArgumentError rref(A2)
    end
end
