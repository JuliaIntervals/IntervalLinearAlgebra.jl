@testset "Verified linear solver" begin

    n = 4
    Arat = reshape(1 .//(1:n^2), n, n)
    brat = Arat*fill(1//1, n)
    Afloat = float(Arat)
    bfloat = float(brat)

    x, cert = epsilon_inflation(Afloat, bfloat)

    @test all(ones(n) .∈ x)
    @test cert

    Ain = convert.(IA.Interval{Float64}, IA.Interval.(Arat, Arat))
    bin = convert.(IA.Interval{Float64}, IA.Interval.(brat, brat))

    x, cert = epsilon_inflation(Ain, bin)

    @test all(ones(n) .∈ x)
    @test cert

    # big float test
    Abig = BigFloat.(Arat)
    bbig = BigFloat.(brat)

    x, cert = epsilon_inflation(Abig, bbig)

    @test cert
    @test all(diam.(x) .< 1e-50)
    @test all(ones(n) .∈ x)

    # case when should not be possible to certify
    A = [1..2 1..4;0..1 0..1]
    b = A*[1; 1]

    _, cert = epsilon_inflation(A, b)
    @test !cert
end
