# Routines to classify interval matrices

"""
    is_strongly_regular(A::AbstractMatrix{T}) where {T<:Interval}

Tests whether an interval matrix ``A`` is strongly regular, i.e. if ``A_c^{-1}A`` is
an H-matrix. For more details see section 4.6 of [[HOR19]](@ref).

### Examples
```jldoctest
julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> is_strongly_regular(A)
true

julia> A = [0..2 1..1;-1.. -1 0..2]
2×2 Matrix{Interval{Float64}}:
   [0, 2]  [1, 1]
 [-1, -1]  [0, 2]

julia> is_strongly_regular(A)
false
```
"""
function is_strongly_regular(A)
    m, n = size(A)
    m == n || return false

    Ac = mid.(A)
    rank(Ac) == n || return false
    return is_H_matrix(inv(Ac)*A)
end

"""
    is_H_matrix(A::AbstractMatrix{T}) where {T<:Interval}

Tests whether an interval matrix A is an H-matrix, by testing that ``⟨A⟩^{-1}e>0``,
where ``e=[1, 1, …, 1]ᵀ``. Note that in practice it tests that a
_floating point approximation_ of ``⟨A⟩^{-1}e`` satisfies the condition.
For more details see section 4.4 of [[HOR19]](@ref).

### Examples
```jldoctest
julia> A = [2..4 -1..1; -1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> is_H_matrix(A)
true

julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> is_H_matrix(A)
false
```
"""
function is_H_matrix(A)
    m, n = size(A)
    m == n || return false

    compA = comparison_matrix(A)
    rank(compA) == n || return false
    return all(compA\ones(n) .> 0)
end


"""
    is_strictly_diagonally_dominant(A::AbstractMatrix{T}) where {T<:Interval}

Checks whether an interval matrix ``A`` is stictly diagonally dominant, that is
if ``mig(Aᵢᵢ) > ∑_{k ≠ i} mag(Aᵢₖ)`` for ``i=1,…,n``.
For more details see section 4.5 of [[HOR19]](@ref).

### Examples
```jldoctest
julia> A = [2..4 -1..1; -1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> is_strictly_diagonally_dominant(A)
true

julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> is_strictly_diagonally_dominant(A)
false
```
"""
function is_strictly_diagonally_dominant(A)
    m, n = size(A)
    m == n || return false

    @inbounds for i=1:m
        sum_mag = sum(Interval(mag(A[i, k])) for k=1:n if k ≠ i)
        mig(A[i, i]) ≤ inf(sum_mag) && return false

    end
      return true
end

"""
    is_Z_matrix(A::AbstractMatrix{T}) where {T<:Interval}

Checks whether the interval matrix ``A`` is a Z-matrix, that is whether ``Aᵢⱼ≤0`` for all
``i≠j``. For more details see section 4.2 of [[HOR19]](@ref).

### Examples
```jldoctest
julia> A = [2..4 -2.. -1; -2.. -1 2..4]
2×2 Matrix{Interval{Float64}}:
   [2, 4]  [-2, -1]
 [-2, -1]    [2, 4]

julia> is_Z_matrix(A)
true

julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> is_Z_matrix(A)
false
```
"""
function is_Z_matrix(A)

    m, n = size(A)
    m == n || return false

    @inbounds for j in 1:n
        for i in 1:m
            if i != j && sup(A[i, j]) > 0
                return false
            end
        end
    end
    return true
end

"""
    is_M_matrix(A::AbstractMatrix{T}) where {T}

Checks whether the interval matrix ``A`` is an M-matrix, that is a Z-matrix with
non-negative inverse. For more details see section 4.2 of [[HOR19]](@ref).

### Examples
```jldoctest
julia> A = [2..2 -1..0; -1..0 2..2]
2×2 Matrix{Interval{Float64}}:
  [2, 2]  [-1, 0]
 [-1, 0]   [2, 2]

julia> is_M_matrix(A)
true

julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> is_M_matrix(A)
false
```
"""
function is_M_matrix(A)

    is_Z_matrix(A) || return false

    Ainf = inf.(A)
    e = ones(size(A, 1))
    u = Ainf\e
    all(u .> 0) || return false
    return all(A*u .> 0)
end
