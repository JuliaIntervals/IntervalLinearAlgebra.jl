@testset "Matrix multiplication" begin

    # test default settings
    @test get_multiplication_mode() == Dict(:multiplication => :fast)

    A = [2..4 -2..1; -1..2 2..4]
    imA = im*A

    set_multiplication_mode(:slow)
    @test all(isequal_interval.(A * A, [0..18 -16..8; -8..16 0..18]))
    @test all(isequal_interval.(imA * imA, -1*[0..18 -16..8; -8..16 0..18]))

    set_multiplication_mode(:rank1)
    @test all(isequal_interval.(A * A, [0..18 -16..8; -8..16 0..18]))
    @test all(isequal_interval.(imA * imA, -1*[0..18 -16..8; -8..16 0..18]))

    set_multiplication_mode(:fast)
    @test all(isequal_interval.(A * A, [-2..19.5 -16..10; -10..16 -2..19.5]))
    @test all(isequal_interval.(A * mid.(A), [5..12.5 -8..2; -2..8 5..12.5]))
    @test all(isequal_interval.(mid.(A) * A, [5..12.5 -8..2; -2..8 5..12.5]))

    @test all(isequal_interval.(imA * imA, -1*[-2..19.5 -16..10; -10..16 -2..19.5]))
    @test all(isequal_interval.(mid.(A) * imA, im*[5..12.5 -8..2; -2..8 5..12.5]))
    @test all(isequal_interval.(imA * mid.(A), im*[5..12.5 -8..2; -2..8 5..12.5]))
end
