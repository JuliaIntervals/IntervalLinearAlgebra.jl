struct AffineParametricArray{T, N, MT<:AbstractArray{T, N}, VAR} <: AbstractArray{T, N}
    coeffs::Vector{MT}
    vars::Vector{VAR}
end

const AffineParametricMatrix{T, MT} = AffineParametricArray{T, 2, MT} where {T, MT<:AbstractMatrix{T}}

const AffineParametricVector{T, VT} = AffineParametricArray{T, 1, VT} where {T, VT <: AbstractVector{T}}

# apas are callable
function (apa::AffineParametricArray)(p::AbstractVector)
    length(p) + 1 == length(apa.coeffs) || throw(ArgumentError("dimension mismatch"))
    return sum(apa.coeffs[i] * p[i] for i in eachindex(p)) + apa.coeffs[end]
end

function (apa::AffineParametricArray)()
   return apa.coeffs[end] + sum(apa.coeffs[i] * apa.vars[i] for i in eachindex(apa.vars))
end

# Array interface
IndexStyle(::Type{<:AffineParametricArray}) = IndexLinear()
size(apa::AffineParametricArray) = size(apa.coeffs[1])

getindex(apa::AffineParametricArray, idx...) = getindex(apa(), idx...)

# operations
function ==(apa1::AffineParametricArray, apa2::AffineParametricArray)
    return apa1.vars == apa2.vars && apa1.coeffs == apa2.coeffs
end

# unary operations
+(apa::AffineParametricArray) = apa
-(apa::AffineParametricArray) = AffineParametricArray(-apa.coeffs, apa.vars)

# addition subtraction
for op in (:+, :-)
    @eval function $op(apa1::AffineParametricArray, apa2::AffineParametricArray)
        return AffineParametricArray($op(apa1.coeffs, apa2.coeffs), apa1.vars)
    end

    @eval function $op(apa::AffineParametricArray{T, N, MT, VAR}, B::MS) where {T, N, MT, VAR, S, MS<:AbstractArray{S, N}}
        MTS = promote_type(MT, MS)
        coeffs = similar(apa.coeffs, MTS)
        coeffs .= apa.coeffs
        coeffs[end] = $op(coeffs[end], B)
        return AffineParametricArray(coeffs, apa.vars)
    end

    @eval function $op(B::MS, apa::AffineParametricArray{T, N, MT, VAR}) where {T, N, MT, VAR, S, MS<:AbstractArray{S, N}}
        MTS = promote_type(MT, MS)
        coeffs = similar(apa.coeffs, MTS)
        coeffs .= $op(apa.coeffs)
        coeffs[end] = $op(B, apa.coeffs[end])
        return AffineParametricArray(coeffs, apa.vars)
    end
end

# multiplication, backslash
for op in (:*, :\)
    @eval function $op(A::AbstractMatrix, apa::AffineParametricArray)
        coeffs = [$op(A, coeff) for coeff in apa.coeffs]
        return AffineParametricArray(coeffs, apa.vars)
    end
end

function *(apa::AffineParametricMatrix, B::AbstractMatrix)
    coeffs = [coeff*B for coeff in apa.coeffs]
    return AffineParametricArray(coeffs, apa.vars)
end

function *(apa::AffineParametricMatrix, v::AbstractVector)
    coeffs = [coeff*v for coeff in apa.coeffs]
    return AffineParametricArray(coeffs, apa.vars)
end

*(apa::AffineParametricArray, n::Number) = AffineParametricArray(apa.coeffs * n, apa.vars)
*(n::Number, apa::AffineParametricArray) = AffineParametricArray(apa.coeffs * n, apa.vars)
/(apa::AffineParametricArray, n::Number) = AffineParametricArray(apa.coeffs / n, apa.vars)

# with AffineExpression
function AffineParametricArray(A::AbstractArray{<:AffineExpression}, vars::Vector{<:AffineExpression})
    ncoeffs = length(vars) + 1
    coeffs = [map(a -> a.coeffs[i], A) for i in 1:ncoeffs]
    return AffineParametricArray(coeffs, vars)
end

function AffineParametricArray(A::AbstractArray{<:Number}, vars::Vector{<:AffineExpression})
    ncoeffs = length(vars) + 1
    coeffs = [zero(A) for _ in 1:ncoeffs]
    coeffs[end] = A
    # vars = [AffineExpression(Vector(c)) for c in eachcol(Matrix{T}(I, ncoeffs, ncoeffs-1))]
    return AffineParametricArray(coeffs, vars)
end

function setindex!(apa::AffineParametricArray, ae::AffineExpression, idx...)
    @inbounds for i in eachindex(ae.coeffs)
        setindex!(apa.coeffs[i], ae.coeffs[i], idx...)
    end
end
