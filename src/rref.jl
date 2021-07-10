"""
    rref(A::AbstractMatrix{T}) where {T<:Interval}

Computes the reduced row echelon form of the interval matrix `A` using maximum
mignitude as pivoting strategy.

### Examples

```jldoctest
julia> A = [2..4 -1..1; -1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> rref(A)
2×2 Matrix{Interval{Float64}}:
 [2, 4]  [-1, 1]
 [0, 0]       [1.5, 4.5]
```
"""
function rref(A::AbstractMatrix{T}) where {T<:Interval}
    A1 = copy(A)
    return rref!(A1)
end

"""
    rref!(A::AbstractMatrix{T}) where {T<:Interval}

In-place version of [`rref`](@ref).
"""
function rref!(A::AbstractMatrix{T}) where {T<:Interval}
    m, n = size(A)
    minmn = min(m,n)
    @inbounds for k = 1:minmn

        if k < m
            # find maximum index
            migmax, kp = _findmax_mig(view(A, k:m, k))
            iszero(migmax) && throw(ArgumentError("Could not find a pivot with non-zero mignitude in column $k."))
            kp += k - 1

            # Swap rows k and kp if needed
            k != kp && _swap!(A, k, kp)
        end

        # Scale first column
        _scale!(A, k)

        # Update the rest
        _eliminate!(A, k)
    end
    return A
end


@inline function _findmax_mig(v)

    @inbounds begin
        migmax = mig(first(v))
        kp = firstindex(v)

        for (i, vi) in enumerate(v)
            migi = mig(vi)
            if migi > migmax
                kp = i
                migmax = migi
            end
        end
    end
    return migmax, kp
end


@inline function _swap!(A, k, kp)
    @inbounds for i = 1:size(A, 2)
        tmp = A[k,i]
        A[k,i] = A[kp,i]
        A[kp,i] = tmp
    end
end

@inline function _scale!(A, k)
    @inbounds begin
        Akkinv = inv(A[k,k])
        for i = k+1:size(A, 1)
            A[i,k] *= Akkinv
        end
    end
end

@inline function _eliminate!(A, k)
    m, n = size(A)
    @inbounds begin
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
