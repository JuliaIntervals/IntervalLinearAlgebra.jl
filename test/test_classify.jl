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
