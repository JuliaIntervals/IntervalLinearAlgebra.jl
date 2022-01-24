"""
    AffineParametricArray{T, N, MT<:AbstractArray{T, N}}

Array whose elements have an affine dependency on a set of parameteres ``p₁, p₂, …, pₙ``.

### Fields

- coeffs::Vector{MT} -- vector of arrays, corresponds to the coefficients of each variable.

### Example

```jldoctest
julia> @affinevars x y z
3-element Vector{AffineExpression{Int64}}:
 x
 y
 z

julia> A = AffineParametricArray([x+y x-1;x+y+z 1])
2×2 AffineParametricMatrix{Int64, Matrix{Int64}}:
 x+y    x-1
 x+y+z  1
```
"""
struct AffineParametricArray{T, N, MT<:AbstractArray{T, N}} <: AbstractArray{T, N}
    coeffs::Vector{MT}
end

const AffineParametricMatrix{T, MT} = AffineParametricArray{T, 2, MT} where {T, MT<:AbstractMatrix{T}}

const AffineParametricVector{T, VT} = AffineParametricArray{T, 1, VT} where {T, VT <: AbstractVector{T}}

# apas are callable
function (apa::AffineParametricArray)(p)
    length(p) + 1 == length(apa.coeffs) || throw(ArgumentError("dimension mismatch"))
    return sum(apa.coeffs[i] * p[i] for i in eachindex(p)) + apa.coeffs[end]
end

# Array interface
IndexStyle(::Type{<:AffineParametricArray}) = IndexLinear()
size(apa::AffineParametricArray) = size(apa.coeffs[1])

function getindex(apa::AffineParametricArray, idx...)
    nvars = length(_vars_dict[:vars])
    vars = [AffineExpression(Vector(c)) for c in eachcol(Matrix(I, nvars + 1, nvars))]
    tmp =  sum(getindex(c, idx...) * v for (v, c) in zip(vars, apa.coeffs[1:end-1]))
    return getindex(apa.coeffs[end], idx...) + tmp
end

function setindex!(apa::AffineParametricArray, ae::AffineExpression, idx...)
    @inbounds for i in eachindex(ae.coeffs)
        setindex!(apa.coeffs[i], ae.coeffs[i], idx...)
    end
end

function setindex!(apa::AffineParametricArray, num::Number, idx...)
    setindex!(apa.coeffs[end], num, idx...)
end

==(apa1::AffineParametricArray, apa2::AffineParametricArray) = apa1.coeffs == apa2.coeffs

# unary operations
+(apa::AffineParametricArray) = apa
-(apa::AffineParametricArray) = AffineParametricArray(-apa.coeffs)

# addition subtraction
for op in (:+, :-)
    @eval function $op(apa1::AffineParametricArray, apa2::AffineParametricArray)
        return AffineParametricArray($op(apa1.coeffs, apa2.coeffs))
    end

    @eval function $op(apa::AffineParametricArray{T, N, MT}, B::MS) where {T, N, MT, S, MS<:AbstractArray{S, N}}
        MTS = promote_type(MT, MS)
        coeffs = similar(apa.coeffs, MTS)
        coeffs .= apa.coeffs
        coeffs[end] = $op(coeffs[end], B)
        return AffineParametricArray(coeffs)
    end

    @eval function $op(B::MS, apa::AffineParametricArray{T, N, MT}) where {T, N, MT, S, MS<:AbstractArray{S, N}}
        MTS = promote_type(MT, MS)
        coeffs = similar(apa.coeffs, MTS)
        coeffs .= $op(apa.coeffs)
        coeffs[end] = $op(B, apa.coeffs[end])
        return AffineParametricArray(coeffs)
    end
end

# multiplication, backslash
for op in (:*, :\)
    @eval function $op(A::AbstractMatrix, apa::AffineParametricArray)
        coeffs = [$op(A, coeff) for coeff in apa.coeffs]
        return AffineParametricArray(coeffs)
    end
end

function *(apa::AffineParametricMatrix, B::AbstractMatrix)
    coeffs = [coeff*B for coeff in apa.coeffs]
    return AffineParametricArray(coeffs)
end

function *(apa::AffineParametricMatrix, v::AbstractVector)
    coeffs = [coeff*v for coeff in apa.coeffs]
    return AffineParametricArray(coeffs)
end

*(apa::AffineParametricArray, n::Number) = AffineParametricArray(apa.coeffs * n)
*(n::Number, apa::AffineParametricArray) = AffineParametricArray(apa.coeffs * n)
/(apa::AffineParametricArray, n::Number) = AffineParametricArray(apa.coeffs / n)

# with AffineExpression
function AffineParametricArray(A::AbstractArray{<:AffineExpression})
    ncoeffs = length(first(A).coeffs)
    coeffs = [map(a -> a.coeffs[i], A) for i in 1:ncoeffs]
    return AffineParametricArray(coeffs)
end

function AffineParametricArray(A::AbstractArray{<:Number})
    ncoeffs = length(_vars_dict[:vars]) + 1
    coeffs = [zero(A) for _ in 1:ncoeffs]
    coeffs[end] = A
    return AffineParametricArray(coeffs)
end
