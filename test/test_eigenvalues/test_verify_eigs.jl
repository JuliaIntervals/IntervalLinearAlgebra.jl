@testset "verify perron frobenius " begin
    A = [1 2;3 4]
    ρ = bound_perron_frobenius_eigenvalue(A)
    @test  (5+sqrt(big(5)))/2 ≤ ρ
end

@testset "verified eigenvalues" begin
    n = 5 # matrix size

    # symmetric case
    ev = sort(randn(n))
    D = Diagonal(ev)
    Q, _ = qr(rand(n, n))
    A = Symmetric(IA.Interval.(Q) * D * IA.Interval.(Q)')

    evals, evecs, cert = verify_eigen(A)
    @test all(cert)
    @test all(ev .∈ evals)


    # real eigenvalues case
    P = rand(n, n)
    Pinv, _ = epsilon_inflation(P, Diagonal(ones(n)))
    A = IA.Interval.(P) * IA.Interval.(D) * Pinv

    evals, evecs, cert = verify_eigen(A)
    @test all(cert)
    @test all(ev .∈ evals)

    # test complex eigenvalues
    ev = sort(rand(Complex{Float64}, n), by = x -> (real(x), imag(x)))
    A = IA.Interval.(P) * Matrix(Diagonal(ev)) * Pinv

    evals, evecs, cert = verify_eigen(A)
    @test all(cert)
    @test all(ev .∈ evals)
end
