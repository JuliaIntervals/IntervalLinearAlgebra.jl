@testset "Matrix multiplication" begin

    # test default settings
    @test get_multiplication_mode() == Dict(:multiplication => :fast)

    A = [2..4 -2..1; -1..2 2..4]
    set_multiplication_mode(:slow)
    @test A*A == [0..18 -16..8; -8..16 0..18]

    set_multiplication_mode(:rank1)
    @test A*A == [0..18 -16..8; -8..16 0..18]

    set_multiplication_mode(:fast)
    @test A*A == [-2..19.5 -16..10; -10..16 -2..19.5]
end
