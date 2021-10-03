@testset "linear expressions" begin
    @linvars x y z

    @test x isa IntervalLinearAlgebra.AffineExpression{Int}
    @test x.coeffs == [1, 0, 0, 0]
    @test y.coeffs == [0, 1, 0, 0]
    @test z.coeffs == [0, 0, 1, 0]

    p1 = x - y + 3z
    @test IntervalLinearAlgebra._tostring(p1) == "x-y+3z"
    p2 = y + x - z - 2
    @test IntervalLinearAlgebra._tostring(p2) == "x+y-z-2"

    @test +p1 == p1
    @test -p1 == -x + y - 3z

    psum = p1 + p2
    pdiff = p1 - p2
    psumnumleft = 1 + p1
    psumnumright = p1 + 1
    pdiffnumright = p1 - 1
    pdiffnumleft = 1 - p1
    pprodleft = 2 * p1
    pprodright = p1 * 2
    pdiv = p1 / 2

    @test psum.coeffs == [2, 0, 2, -2]
    @test pdiff.coeffs == [0, -2, 4, 2]
    @test psumnumleft.coeffs == [1, -1, 3, 1]
    @test psumnumright.coeffs == [1, -1, 3, 1]
    @test pdiffnumright.coeffs ==  [1, -1, 3, -1]
    @test pdiffnumleft.coeffs == [-1, 1, -3, 1]
    @test pprodleft.coeffs == [2, -2, 6, 0]
    @test pprodright.coeffs == [2, -2, 6, 0]
    @test pdiv.coeffs == [0.5, -0.5, 1.5, 0]

end
