"""
    interval_norm(A::Matrix{Interval})

computes the infinity norm of interval matrix A.
"""
interval_norm(A::AbstractMatrix) = opnorm(mag.(A), Inf) # TODO: proper expand norm function and add 1-norm
interval_norm(v::AbstractVector) = maximum(mag.(v))
# ? use manual loops instead

"""
Preconditions the interval system Ax = b by multipling by the (approximate) inverse of Ac,
that is the midpoint of A.
"""
function precondition(A, b)
    Ac_inv = inv(mid.(A))
    return Ac_inv*A, Ac_inv*b
end


"""
    enclose(A::Matrix{Interval}, b::Vector{Interval})

Computes an initial enclosure Σ so that x ⊆ Σ, where x is the solution of the interval
system Ax = b. It assumes the system has already been preconditioned or does not require
preconditioning. See proposition 5.14 of [1] (page 51)
"""
function enclose(A::StaticMatrix{N, N, T}, b::StaticVector{N, T}) where {N, T}
    A1 = Diagonal(ones(N)) - A
    e = interval_norm(b)/(1 - interval_norm(A1))
    x0 = MVector{N, T}(fill(-e..e, N))
    return x0
end

function enclose(A, b)
    n = length(b)
    A1 = Diagonal(ones(n)) - A
    e = interval_norm(b)/(1 - interval_norm(A1))
    x0 = fill(-e..e, n)
    return x0
end

"""
    comparison_matrix(A::Matrix{Interval})

Computes the comparison matrix ⟨A⟩ of the given matrix A according to the definition
⟨A⟩_ii = mig(A_ii)
⟨A⟩_ij = -mag(A_ij)
See definition 4.16 of [1] (page 33)
"""
function comparison_matrix(A::SMatrix)
    n = size(A, 1)
    compA = -mag.(A)
    @inbounds for (i, idx) in enumerate(diagind(A))
        compA = setindex(compA, mig(A[i, i]), idx)
    end
    return compA
end

function comparison_matrix(A)
    n = size(A, 1)
    compA = -mag.(A)
    @inbounds for i in 1:n
        compA[i, i] = mig(A[i, i])
    end
    return compA
end
