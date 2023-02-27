module NumericalTest

export rounding_test

function test_matrix(k)
    A = zeros(Float64, (k, k))
    A[:, end] = fill(2^(-53), k)
    for i in 1:k-1
        A[i,i] = 1.0
    end
    return A
end

using LinearAlgebra
"""
    rounding_test(n, k)

Let `u=fill(2^(-53), k-1)` and let A be the matrix
[I u;
0 2^(-53)]

This test checks the result of A*A' in different rounding modes,
running BLAS on `n` threads
"""
function rounding_test(n,k)
    
    
    BLAS.set_num_threads( n )
    A = test_matrix( k )
    B = setrounding(Float64, RoundUp) do
        BLAS.gemm('N', 'T', 1.0, A, A)
    end

    return all([B[i,i]==nextfloat(1.0) for i in 1:k-1])
end


end