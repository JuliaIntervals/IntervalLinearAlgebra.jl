@testset "verify perron frobenius " begin
    A = [1 2;3 4]
    ρ = bound_perron_frobenius_eigenvalue(A)
    @test  (5+sqrt(big(5)))/2 ≤ ρ
end

@testset "verified eigenvalues of symmetric matrix" begin
    n = 10 # matrix size
    ev = sort(randn(n))
    D = Diagonal(ev)
    Q, _ = qr(rand(n, n))
    A = Symmetric(IA.Interval.(Q) * D * IA.Interval.(Q)')

    F, cert = verify_eigenvalues(A)
    @test all(cert)
    @test all(ev .∈ F.values)
end
