using IntervalLinearAlgebra, StaticArrays, IntervalArithmetic
using Test

int_approx(a, b, tol) =  abs(a.lo - b.lo) < tol  && abs(a.hi - b.hi) < tol

@testset "IntervalLinearAlgebra.jl" begin
    
    @testset "Test linear solvers" begin
        As = @SMatrix [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]
        bs = @SVector [-2..4, 1..8, -4..10, 2..12]

        Am = Matrix(As)
        bm = Vector(bs)
        jac = Jacobi()
        gs = GaussSeidel()
        hbr = HansenBliekRohn()

        for (A, b) in zip([As, Am], [bs, bm])        
            xgs = solve(A, b, gs)
            xjac = solve(A, b, jac)
            xhbr = solve(A, b, hbr)

            @test all(int_approx.(xgs, [-2.6..3.1, -3.9..1.65, -1.48..2.15, -2.35..0.79], 0.01))
            @test all(int_approx.(xjac, [-2.6..3.1, -3.9..1.65, -1.48..2.15, -2.35..0.79], 0.01))

            @test all(int_approx.(xhbr, [-2.5..3.1, -3.9..1.2, -1.4..2.15, -2.35..0.6], 0.01))
        end
    end
end
