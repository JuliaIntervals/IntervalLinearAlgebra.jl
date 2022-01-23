@testset "Affine Parametric Array construction" begin
    @affinevars x y z

    vars = [x, y, z]

    A = AffineParametricArray([x+1 y+2;x+y+z+1 x-z])
    A1 = AffineParametricArray([1 2;3 4])
    B = AffineParametricArray([x+1 y+2;x+y+z+1 x-z+1])

    @test A.coeffs == [[1 0;1 1], [0 1;1 0], [0 0;1 -1], [1 2;1 0]]
    @test A1.coeffs == [[0 0;0 0], [0 0;0 0], [0 0;0 0], [1 2;3 4]]

    @test A[1, 2] == y + 2
    @test A[:, 1] == [x+1, x+y+z+1]

    A1[1, 2] = x
    @test A1 == AffineParametricArray([1 x;3 4])
    A1[:, 1] = [x+1, z-1]
    @test A1 == AffineParametricArray([x+1 x;z-1 4])

    @test A([1..2, 2..3, 3..4]) == [2..3 4..5; 7..10 -3.. -1]
end

@testset "Affine parametric array operations" begin
    @affinevars x y z
    vars = [x, y, z]

    A = AffineParametricArray([x+1 y+2;x+y+z+1 x-z])
    B = AffineParametricArray([x+1 y+2;x+y+z+1 x-z+1])

    @test A == A
    @test A != B

    @test +A == A
    @test -A == AffineParametricArray([-x-1 -y-2; -x-y-z-1 -x+z])
    @test A + B == AffineParametricArray([2x+2 2y+4; 2x+2y+2z+2 2x-2z+1])
    @test A - B == AffineParametricArray([0 0;0 -1])

    C = [2 1;1 1]
    @test A + C == AffineParametricArray([x+3 y+3;x+y+z+2 x-z+1])
    @test C + A == AffineParametricArray([x+3 y+3;x+y+z+2 x-z+1])

    @test A * C == AffineParametricArray([2x+y+4 x+y+3;3x+2y+z+2 2x+y+1])
    @test C * A == AffineParametricArray([3x+y+z+3 x+2y-z+4;2x+y+z+2 x+y-z+2])
    @test C \ A == inv(C) * A

    @test A * 2 == AffineParametricArray([2x+2 2y+4;2x+2y+2z+2 2x-2z])
    @test 2 * A == AffineParametricArray([2x+2 2y+4;2x+2y+2z+2 2x-2z])
    @test A / 2 == AffineParametricArray([0.5x+0.5 0.5y+1;0.5x+0.5y+0.5z+0.5 0.5x-0.5z])
end
