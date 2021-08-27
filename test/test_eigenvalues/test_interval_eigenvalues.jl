@testset "Eigenvalues of interval matrices" begin

    # symmetrix matrix
    A = Symmetric([-1 0 -1..1;
         0 -1 -1..1;
         -1..1 -1..1 0.1])

    ev = eigenbox(A)
    @test interval_isapprox(ev, -2.4143..1.5143; atol=1e-3)

    # real matrix
    A = [-3.. -2 4..5 4..6 -1..1.5;
         -4.. -3 -4.. -3 -4.. -3 1..2;
         -5.. -4 2..3 -5.. -4 -1..0;
         -1..0.1 0..1 1..2 -4..2.5]

    ev = eigenbox(A)
    @test interval_isapprox(real(ev), -8.8221..3.4408; atol=1e-3)
    @test interval_isapprox(imag(ev), -10.7497..10.7497; atol=1e-3)


    # hermitian matrix
    A = Hermitian([1..2 (5..9)+(2..5)*im (3..5)+(2..4)im;
        (5..9)+(-5.. -2)*im 2..3 (7..8)+(6..10)im;
        (3..5)+(-4.. -2)*im (7..8)+(-10.. -6)*im 3..4])

    ev = eigenbox(A)
    @test interval_isapprox(ev, -15.4447..24.3359; atol=1e-3)

    # complex matrix
    A = [(1..2)+(3..4)*im 3..4;1+(2..3)*im 4..5]

    ev = eigenbox(A)
    @test interval_isapprox(real(ev), -1.28812..7.28812; atol=1e-3)
    @test interval_isapprox(imag(ev), -2.04649..5.54649; atol=1e-3)
end
