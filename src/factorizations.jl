"""
    rref(A)

computes the reduced row echelon form of the interval matrix A using maximum mignitude as pivoting strategy.

### Parameters:
A : interval matrix
"""
function rref(A)
    A1 = copy(A)
    return rref!(A1)
end

function rref!(A)
    m, n = size(A)
    minmn = min(m,n)
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if k < m
                migmax = mig(A[k, k])
                for i = k+1:m
                    migi = mig(A[i,k])
                    if migi > migmax
                        kp = i
                        migmax = migi
                    end
                end
                iszero(migmax) && throw(ArgumentError("Could not find a pivot with non-zero mignitude in column $k."))
            end

            if k != kp
                # Interchange
                for i = 1:n
                    tmp = A[k,i]
                    A[k,i] = A[kp,i]
                    A[kp,i] = tmp
                end
            end
            # Scale first column
            Akkinv = inv(A[k,k])
            for i = k+1:m
                A[i,k] *= Akkinv
            end

            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end

            for i = k+1:m
                A[i, k] = zero(eltype(A))
            end
        end
    end
    return A
end
