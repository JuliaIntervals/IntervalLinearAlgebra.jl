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
  [2.0, 4.0]_com  [-2.0, 1.0]_com
 [-1.0, 2.0]_com   [2.0, 4.0]_com

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2.0, 2.0]_com
 [-2.0, 2.0]_com

julia> np = NoPrecondition()
NoPrecondition()

julia> np(A, b)
(Interval{Float64}[[2.0, 4.0]_com [-2.0, 1.0]_com; [-1.0, 2.0]_com [2.0, 4.0]_com], Interval{Float64}[[-2.0, 2.0]_com, [-2.0, 2.0]_com])
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
  [2.0, 4.0]_com  [-2.0, 1.0]_com
 [-1.0, 2.0]_com   [2.0, 4.0]_com

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2.0, 2.0]_com
 [-2.0, 2.0]_com

julia> imp = InverseMidpoint()
InverseMidpoint()

julia> imp(A, b)
(Interval{Float64}[[0.594595, 1.40541]_com [-0.540541, 0.540541]_com; [-0.540541, 0.540541]_com [0.594595, 1.40541]_com], Interval{Float64}[[-0.756757, 0.756757]_com_NG, [-0.756757, 0.756757]_com_NG])
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
  [2.0, 4.0]_com  [-2.0, 1.0]_com
 [-1.0, 2.0]_com   [2.0, 4.0]_com

julia> b = [-2..2, -2..2]
2-element Vector{Interval{Float64}}:
 [-2.0, 2.0]_com
 [-2.0, 2.0]_com

julia> idmp = InverseDiagonalMidpoint()
InverseDiagonalMidpoint()

julia> idmp(A, b)
(Interval{Float64}[[0.666667, 1.33333]_com [-0.666667, 0.333333]_com; [-0.333333, 0.666667]_com [0.666667, 1.33333]_com], Interval{Float64}[[-0.666667, 0.666667]_com_NG, [-0.666667, 0.666667]_com_NG])
```
"""
struct InverseDiagonalMidpoint <: AbstractPrecondition end

function (idmp::InverseDiagonalMidpoint)(A::AbstractMatrix{T},
                                         b::AbstractVector{T}) where {T<:Interval}
    R = inv(Diagonal(mid.(A)))
    return R*A, R*b
end
