@testset "test skalna06" begin

    E = 2e11
    σ = 0.005

    s12 = s21 = s34 = s43 = s45 = s54 = E*σ/sqrt(2)
    s13 = s31 = s24 = s42 = s35 = s53 = E*σ/2

    # dummy variable
    @affinevars s23

    K = AffineParametricArray([s12/2+s13 -s12/2 -s12/2 -s13 0 0 0;
        -s21/2 (s21+s23)/2+s24 (s21-s23)/2 -s23/2 s23/2 -s24 0;
        -s21/2 (s21-s23)/2 (s21+s23)/2 s23/2 -s23/2 0 0;
        -s31 -s23/2 s23/2 s31+(s23+s34)/2+s35 (s34-s23)/2 -s34/2 -s34/2;
        0 s23/2 -s23/2 (s34 - s23)/2 (s34+s23)/2 -s34/2 -s34/2;
        0 -s42 0 -s43/2 -s43/2 s42+(s43+s45)/2 0;
        0 0 0 -s43/2 -s43/2 0 (s43+s45)/2])

    q = [0, 0, -10.0^4, 0, 0, 0, 0]

    s = E*σ/sqrt(2) ± 0.1 * E*σ/sqrt(2)

    _x = [-20, -2.7.. -2.3, -38.91.. -38.52, -5, -34.53.. -33.75, -12.7.. -12.3, -19.77.. -19.37]

    x = solve(K, q, s)/1e-6
    @test all(interval_isapprox(xi, _xi; atol=1e-2) for (xi, _xi) in zip(x, _x))

end
