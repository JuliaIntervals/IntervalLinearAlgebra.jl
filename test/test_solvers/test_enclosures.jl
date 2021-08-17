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

    # test exceptions
    @test_throws DimensionMismatch solve(Am, bm[1:end-1])
    @test_throws DimensionMismatch solve(Am[:, 1:end-1], bm, hbr)
    @test_throws DimensionMismatch solve(Am, bm, gs, NoPrecondition(), [1..2, 3..4])
end

@testset "Reduced Row Echelon Form" begin
    A1 = [1..2 1..2;2..2 3..3]
    @test rref(A1) == [2..2 3..3; 0..0 -2..0.5]

    A2 = fill(0..0, 2, 2)
    @test_throws ArgumentError rref(A2)
end
