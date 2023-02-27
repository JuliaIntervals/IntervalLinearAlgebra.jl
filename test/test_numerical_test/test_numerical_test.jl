@testset "Test Numerical Test" begin
    A = IntervalLinearAlgebra.NumericalTest.test_matrix(4)
    @test A == [1.0 0 0 2^(-53);
                0 1.0 0 2^(-53);
                0 0 1.0 2^(-53);
                0 0 0 2^(-53)]

    # we test the singlethread version, so that CI points 
    # out if setting rounding modes is broken
    @test IntervalLinearAlgebra.NumericalTest.rounding_test(1, 2) == true
end