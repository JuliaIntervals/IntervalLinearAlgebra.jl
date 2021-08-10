"""
    AbstractPrecondition

Abstract type for preconditioners of interval linear systems.
"""
abstract type AbstractPrecondition end

"""
    NoPrecondition <: AbstractPrecondition

Type of the trivial preconditioner which does nothing.

### Notes

- An object of type `NoPrecondition` is a function with method

        (np::NoPrecondition)(A::AbstractMatrix{T},
                             b::AbstractVector{T}) where {T<:Interval}

### Example

```jldoctest
julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]

julia> np = NoPrecondition()
NoPrecondition()

julia> np(A, b)
(Interval{Float64}[[2, 4] [-2, 1]; [-1, 2] [2, 4]], Interval{Float64}[[-2, 2], [-2, 2]])
```
"""
struct NoPrecondition <: AbstractPrecondition end

(np::NoPrecondition)(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T<:Interval} = A, b

"""
    InverseMidpoint <: AbstractPrecondition

Preconditioner that preconditions the linear system ``Ax=b`` with ``A_c^{-1}``,
where ``A_c`` is the midpoint matrix of ``A``.

### Notes

- An object of type `InverseMidpoint` is a function with method

        (imp::InverseMidpoint)(A::AbstractMatrix{T},
                               b::AbstractVector{T}) where {T<:Interval}

### Examples

```jldoctest
julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]

julia> imp = InverseMidpoint()
InverseMidpoint()

julia> imp(A, b)
(Interval{Float64}[[0.594594, 1.40541] [-0.540541, 0.540541]; [-0.540541, 0.540541] [0.594594, 1.40541]], Interval{Float64}[[-0.756757, 0.756757], [-0.756757, 0.756757]])
```
"""
struct InverseMidpoint <: AbstractPrecondition end

function (imp::InverseMidpoint)(A::AbstractMatrix{T},
                                b::AbstractVector{T}) where {T<:Interval}
    R = inv(mid.(A))
    return R*A, R*b
end

"""
    InverseDiagonalMidpoint <: AbstractPrecondition

Preconditioner that preconditions the linear system ``Ax=b`` with the diagonal matrix of
``A_c^{-1}``, where ``A_c`` is the midpoint matrix of ``A``.

### Notes

- An object of type `InverseDiagonalMidpoint` is a function with method

        (idmp::InverseDiagonalMidpoint)(A::AbstractMatrix{T},
                                        b::AbstractVector{T}) where {T<:Interval}

### Example

```jldoctest
julia> A = [2..4 -2..1; -1..2 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-2, 1]
 [-1, 2]   [2, 4]

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]

julia> idmp = InverseDiagonalMidpoint()
InverseDiagonalMidpoint()

julia> idmp(A, b)
(Interval{Float64}[[0.666666, 1.33334] [-0.666667, 0.333334]; [-0.333334, 0.666667] [0.666666, 1.33334]], Interval{Float64}[[-0.666667, 0.666667], [-0.666667, 0.666667]])
```
"""
struct InverseDiagonalMidpoint <: AbstractPrecondition end

function (idmp::InverseDiagonalMidpoint)(A::AbstractMatrix{T},
                                         b::AbstractVector{T}) where {T<:Interval}
    R = inv(Diagonal(mid.(A)))
    return R*A, R*b
end
