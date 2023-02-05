# mid(x::Real) = x #now exported by IntervalArithmetic
mid(x::Complex) = x

"""
    interval_isapprox(a::Interval, b::Interval; kwargs)

Checks whether the intervals ``a`` and ``b`` are approximate equal, that is both their
lower and upper bound are approximately equal.

### Keywords

Same of `Base.isapprox`

### Example

```jldoctest
julia> a = 1..2
[1, 2]

julia> b = a + 1e-10
[1, 2.00001]

julia> interval_isapprox(a, b)
true

julia> interval_isapprox(a, b; atol=1e-15)
false
```
"""
interval_isapprox(a::Interval, b::Interval; kwargs...) = isapprox(a.lo, b.lo; kwargs...) && isapprox(a.hi, b.hi; kwargs...)

"""
    interval_norm(A::AbstractMatrix{T}) where {T<:Interval}

computes the infinity norm of interval matrix ``A``.

### Examples

```jldoctest
julia> A = [2..4 -1..1; -1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> interval_norm(A)
5.0
```
"""
interval_norm(A::AbstractMatrix{T}) where {T<:Interval} = opnorm(mag.(A), Inf) # TODO: proper expand norm function and add 1-norm

"""
    interval_norm(A::AbstractVector{T}) where {T<:Interval}

computes the infinity norm of interval vector ``v``.

### Examples

```jldoctest
julia> b = [-2..2, -3..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-3, 2]

julia> interval_norm(b)
3.0
```
"""
interval_norm(v::AbstractVector{T}) where {T<:Interval} = maximum(mag.(v))
# ? use manual loops instead

"""
    enclose(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T<:Interval}

Computes an enclosure of the solution of the interval linear system ``Ax=b`` using the
algorithm described in sec. 5.7.1 of [[HOR19]](@ref).
"""
function enclose(A::StaticMatrix{N, N, T}, b::StaticVector{N, T}) where {N, T<:Interval}
    C = inv(mid.(A))
    A1 = Diagonal(ones(N)) - C*A
    e = interval_norm(C*b)/(1 - interval_norm(A1))
    x0 = MVector{N, T}(ntuple(_ -> -e..e, Val(N)))
    return x0
end

function enclose(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T<:Interval}
    n = length(b)
    C = inv(mid.(A))
    A1 = Diagonal(ones(n)) - C*A
    e = interval_norm(C*b)/(1 - interval_norm(A1))
    x0 = fill(-e..e, n)
    return x0
end

"""
    comparison_matrix(A::AbstractMatrix{T}) where {T<:Interval}

Computes the comparison matrix ``⟨A⟩`` of the given interval matrix ``A`` according to the
definition ``⟨A⟩ᵢᵢ = mig(Aᵢᵢ)`` and ``⟨A⟩ᵢⱼ = -mag(Aᵢⱼ)`` if ``i≠j``.

### Examples

```jldoctest
julia> A = [2..4 -1..1; -1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> comparison_matrix(A)
2×2 Matrix{Float64}:
  2.0  -1.0
 -1.0   2.0
```
"""
function comparison_matrix(A::SMatrix{N, N, T, M}) where {N, M, T<:Interval}
    n = size(A, 1)
    compA = -mag.(A)
    @inbounds for (i, idx) in enumerate(diagind(A))
        compA = setindex(compA, mig(A[i, i]), idx)
    end
    return compA
end

function comparison_matrix(A::AbstractMatrix{T}) where {T<:Interval}
    n = size(A, 1)
    compA = -mag.(A)
    @inbounds for i in 1:n
        compA[i, i] = mig(A[i, i])
    end
    return compA
end


"""
    Orthants

Iterator to go through all the ``2ⁿ`` vectors of length ``n`` with elements ``±1``.
This is equivalento to going through the orthants of an ``n``-dimensional euclidean space.

### Fields

n::Int -- dimension of the vector space

### Example

```jldoctest
julia> for or in Orthants(2)
       @show or
       end
or = [1, 1]
or = [-1, 1]
or = [1, -1]
or = [-1, -1]
```
"""
struct Orthants
    n::Int
end

Base.eltype(::Type{Orthants}) = Vector{Int}
Base.length(O::Orthants) = 2^(O.n)

function Base.iterate(O::Orthants, state=1)
    state > 2 ^ O.n && return nothing

    vec = -2*digits(state-1, base=2, pad=O.n) .+ 1
    return (vec, state+1)
end

function Base.getindex(O::Orthants, i::Int)
    1 <= i <= length(O) || throw(BoundsError(O, i))
    return -2*digits(i-1, base=2, pad=O.n) .+ 1
end

Base.firstindex(O::Orthants) = 1
Base.lastindex(O::Orthants) = length(O)


_unchecked_interval(x::Real) = Interval(x)
_unchecked_interval(x::Complex) = Interval(real(x)) + Interval(imag(x)) * im
