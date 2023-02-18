@testset "Test Numerical Test"
    A = IntervalLinearAlgebra.NumericalTest.test_matrix(4)
    @test A == [1.0 0 0 2^(-53);
                0 1.0 0 2^(-53);
                0 0 1.0 2^(-53);
                0 0 0 2^(-53)]

    @test IntervalLinearAlgebra.NumericalTest.rounding_test(4) == true
end