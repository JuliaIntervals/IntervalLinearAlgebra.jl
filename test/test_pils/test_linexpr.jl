@testset "linear expressions variables" begin
    @linvars x y z

    @test x isa AffineExpression{Int}
    @test x.coeffs == [1, 0, 0, 0]
    @test y.coeffs == [0, 1, 0, 0]
    @test z.coeffs == [0, 0, 1, 0]

    @linvars x[1:5]

    @test x1 isa AffineExpression{Int}
    @test x1.coeffs == [1, 0, 0, 0, 0, 0]
    @test x5.coeffs == [0, 0, 0, 0, 1, 0]

    @linvars x

    @test x isa AffineExpression{Int}
    @test x.coeffs == [1, 0]
end

@testset "linear expressions operations" begin
    @linvars x y z

    p1 = x - y + 3z
    @test IntervalLinearAlgebra._tostring(p1) == "x-y+3z"
    p2 = y + x - z - 2
    @test IntervalLinearAlgebra._tostring(p2) == "x+y-z-2"
    @test IntervalLinearAlgebra._tostring(p1 - p1) == "0"

    @test +p1 == p1
    @test -p1 == -x + y - 3z
    @test p2([1, 1, 1]) == -1

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

    @test p2 + 0.5 == x + y - z - 1.5

    A = [1 2;3 4]
    @test x * A == [x 2x;3x 4x]
    @test A * x == [x 2x;3x 4x]

    @test zero(p1) == AffineExpression(0)
    @test one(p1) == AffineExpression(1)
    @test zero(AffineExpression{Int}) == AffineExpression(0)
    @test one(AffineExpression{Int}) == AffineExpression(1)
end

@testset "linear expressions conversions" begin
    @test promote_type(AffineExpression{Int}, Float64) == AffineExpression{Float64}
    @test promote_type(AffineExpression{Int}, AffineExpression{Float64}) == AffineExpression{Float64}

    @linvars x y

    p1 = x + y + 1
    a, b = promote(p1, 1.5)

    @test a isa AffineExpression{Float64}
    @test b isa AffineExpression{Float64}

    @test a == x + y + 1
    @test b == 1.5

    @test convert(AffineExpression, 1.2) isa AffineExpression{Float64}
    @test convert(AffineExpression, 1.2) == 1.2
end
