@testset "precondition" begin
    A = [2..4 -2..1; -1..2 2..4]
    b = [-2..2, -2..2]

    np = NoPrecondition()
    idp = InverseDiagonalMidpoint()
    imp = InverseMidpoint()

    A1, b1 = np(A, b)
    @test A1 == A && b1 == b

    A2, b2 = idp(A, b)
    Acorrect = [2/3..4/3 -2/3..1/3; -1/3..2/3 2/3..4/3]
    bcorrect = [-2/3..2/3, -2/3..2/3]
    @test all(interval_isapprox.(A2, Acorrect)) && all(interval_isapprox.(b2, bcorrect))

    A3, b3 = imp(A, b)
    Acorrect = [22/37..52/37 -20/37..20/37;-20/37..20/37 22/37..52/37]
    bcorrect = [-28/37..28/37, -28/37..28/37]
    @test all(interval_isapprox.(A3, Acorrect)) && all(interval_isapprox.(b3, bcorrect))
end
