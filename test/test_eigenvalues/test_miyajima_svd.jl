@testset "verify singular value perron frobenius " begin
    A = [0.5 0.5;0.3 0.7]
    ρ = IntervalLinearAlgebra._bound_perron_frobenius_singularvalue(A)
    @test  1 ≤ ρ
end

@testset "verified svd" begin
    A = [0.9..1.1 0 0 0 2; 0 0 3 0 0; 0 0 0 0 0; 0 2 0 0 0]
    Σ = IntervalLinearAlgebra.svdbox(A, IntervalLinearAlgebra.M1())

    @test all([3 ∈ Σ[1], 5 ∈ Σ[2]^2, 2 ∈ Σ[3], 0 ∈ Σ[4]])
end
