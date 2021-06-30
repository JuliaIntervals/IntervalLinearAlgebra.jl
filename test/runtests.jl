using IntervalLinearAlgebra, StaticArrays, IntervalConstraintProgramming
using Test

@testset "IntervalLinearAlgebra.jl" begin
    
    @testset "Test linear solvers" begin
        As = @SMatrix [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]
        bs = @SVector [-2..4, 1..8, -4..10, 2..12]

        Am = Matrix(As)
        bm = Vector(bs)
        jac = Jacobi()
        gs = GaussSeidel()
        hbr = HansenBliekRohn()
        kra = Krawczyk()

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
    end
    
    @testset "oettli-präger method" begin
        
        A = [2..4 -2..1; -1..2 2..4]
        b = [-2..2, -2..2]
        vars = (:x, :y)

        X = IntervalBox(-14..14, 2)
        p = oettli(A, b, X, vars)

        for pnt in [[-4, -3], [3, -4], [4, 3], [-3, 4]]
            @test any(pnt ∈ x for x in p.boundary)
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
end

