"""
    interval_norm(A::Matrix{Interval})

computes the infinite norm of interval matrix A
"""
interval_norm(A::AbstractMatrix) = maximum(sum(mag.(A); dims=2)) # TODO: proper expand norm function and add 1-norm
interval_norm(v::AbstractVector) = maximum(mag.(v))
"""
Preconditions the interval system Ax = b by multipling by the (approximate) inverse of Ac,
that is the midpoint of A.
"""
function precondition(A, b)
    Ac = mid.(A)
    return Ac\A, Ac\b
end


"""
    enclose(A::Matrix{Interval}, b::Vector{Interval})

Computes an initial enclosure Σ so that x ⊆ Σ, where x is the solution of the interval
system Ax = b.
"""
function enclose(A::SMatrix{N, N, T, M}, b::SVector{N, T}) where {N, T, M}
    A1 = I - A
    e = interval_norm(b)/(1 - interval_norm(A1))
    x0 = MVector{N, T}(fill(-e..e, N))
    return x0
end
function enclose(A::MMatrix{N, N, T, M}, b::MVector{N, T}) where {N, T, M}
    A1 = I - A
    e = interval_norm(b)/(1 - interval_norm(A1))
    x0 = MVector{N, T}(fill(-e..e, N))
    return x0
end

function enclose(A, b)
    A1 = I - A
    e = interval_norm(b)/(1 - interval_norm(A1))
    x0 = fill(-e..e, length(b)  )
    return x0
end

"""
    comparison_matrix(A::Matrix{Interval})

Computes the comparison matrix ⟨A⟩ of the given matrix A according to the definition
⟨A⟩_ii = mig(A_ii)
⟨A⟩_ij = -mag(A_ij)
"""
function comparison_matrix(A)
    n = size(A, 1)
    compA = -mag.(A)
    @inbounds for (i, idx) in enumerate(diagind(A))
        compA = setindex(compA, mig(A[i, i]), idx)
    end
    return compA
end
