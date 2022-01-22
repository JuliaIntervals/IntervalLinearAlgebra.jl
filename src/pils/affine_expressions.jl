const _vars_dict = Dict(:vars => Symbol[])

"""
    AffineExpression{T}

Data structure to represent affine expressions, such as ``x+2y+z+4``.

### Examples

```jldoctest
julia> @linvars x y z
3-element Vector{AffineExpression{Int64}}:
 x
 y
 z

julia> p1 = x + 2y + z + 4
x+2y+z+4
```
"""
struct AffineExpression{T}
    coeffs::Vector{T}
end

function AffineExpression(x::T) where {T<:Number}
    c = zeros(promote_type(T, Int), length(_vars_dict[:vars]) + 1)
    c[end] = x
    return AffineExpression(c)
end

## VARIABLES CONSTRUCTION

_get_vars(x::Symbol) = [x]

function _get_vars(x...)

    if length(x) == 1
        s = x[1].args[1]
        start = x[1].args[2].args[2]
        stop = x[1].args[2].args[3]
        return [Symbol(s, i) for i in start:stop]
    else
        return collect(x)
    end
end

"""
    @affinevars(x...)

Macro to construct the variables used to represent [`AffineExpression`](@ref).

### Examples

```jldoctest
julia> @affinevars x
1-element Vector{AffineExpression{Int64}}:
 x

julia> @affinevars x y z
3-element Vector{AffineExpression{Int64}}:
 x
 y
 z

julia> @affinevars x[1:4]
4-element Vector{AffineExpression{Int64}}:
 x1
 x2
 x3
 x4
```
"""
macro affinevars(x...)
    vars = _get_vars(x...)
    _vars_dict[:vars] = vars

    ex = quote end
    vars_ex = Expr(:vect)
    for (i, s) in enumerate(vars)
        c = zeros(Int, length(vars) + 1)
        c[i] = 1
        push!(ex.args, :($(esc(s)) = AffineExpression($c)))
        push!(vars_ex.args, s)
    end
    push!(ex.args, esc(vars_ex))
    return ex
end

(ae::AffineExpression)(p::Vector{<:Number}) = dot(ae.coeffs[1:end-1], p) + ae.coeffs[end]

## printing
function _tostring(ae::AffineExpression)
    iszero(ae.coeffs) && return "0"
    str = ""
    @inbounds for (i, x) in enumerate(_vars_dict[:vars])
        c = ae.coeffs[i]
        iszero(c) && continue
        sgn = c > 0 ? "+" : "-"
        c = abs(c) == 1 ? "" : "$(abs(c))"
        str *= sgn * c * "$(x)"
    end

    c = last(ae.coeffs)
    if !iszero(c)
        sgn = c > 0 ? "+" : ""
        str *= sgn * "$(c)"
    end
    str = (str[1] == '+' ? str[2:end] : str)
    return str
end

show(io::IO, ae::AffineExpression) = print(io, _tostring(ae))


#########################
# BASIC FUNCTIONS       #
#########################

function zero(::AffineExpression{T}) where {T}
    return AffineExpression(zeros(T, length(_vars_dict[:vars]) + 1))
end

function zero(::Type{AffineExpression{T}}) where {T}
    return AffineExpression(zeros(T, length(_vars_dict[:vars]) + 1))
end

one(::Type{AffineExpression{T}}) where {T} = zero(AffineExpression{T}) + 1
one(::AffineExpression{T}) where {T} = zero(AffineExpression{T}) + 1

# coefficients(ae::AffineExpression) = ae.coeffs[1:end-1]
# coefficient(ae::AffineExpression, i) = ae.coeffs[i]
# offset(ae::AffineExpression) = ae.coeffs[end]

#########################
# ARITHMETIC OPERATIONS #
#########################

for op in (:+, :-)
    @eval function $op(ae1::AffineExpression, ae2::AffineExpression)
        return AffineExpression($op(ae1.coeffs, ae2.coeffs))
    end

    @eval function $op(ae::AffineExpression{T}, n::S) where {T<:Number, S<:Number}
        TS = promote_type(T, S)
        c = similar(ae.coeffs, TS)
        c .= ae.coeffs
        c[end] = $op(c[end], n)
        return AffineExpression(c)
    end
end

+(ae::AffineExpression) = ae
-(ae::AffineExpression) = AffineExpression(-ae.coeffs)

function Base.:(*)(ae1::AffineExpression, n::Number)
    return AffineExpression(ae1.coeffs * n)
end

function Base.:(/)(ae1::AffineExpression, n::Number)
    return AffineExpression(ae1.coeffs / n)
end

for op in (:+, :*)
    @eval $op(n::Number, ae::AffineExpression) = $op(ae, n)
end

-(n::Number, ae::AffineExpression) = n + (-ae)

==(ae1::AffineExpression, ae2::AffineExpression) = ae1.coeffs == ae2.coeffs
==(ae1::AffineExpression, n::Number) = iszero(ae1.coeffs[1:end-1]) && ae1.coeffs[end] == n

function *(ae::AffineExpression, A::AbstractArray{T, N}) where {T<:Number, N}
    return map(x -> x * ae, A)
end

function *(A::AbstractArray{T, N}, ae::AffineExpression) where {T<:Number, N}
    return map(x -> x * ae, A)
end

## Convetion and promotion

function promote_rule(::Type{AffineExpression{T}}, ::Type{S}) where {T<:Number, S<:Number}
    AffineExpression{promote_type(T, S)}
end

function promote_rule(::Type{AffineExpression{T}}, ::Type{AffineExpression{S}}) where {T<:Number, S<:Number}
    AffineExpression{promote_type(T, S)}
end

convert(::Type{AffineExpression}, x::Number) = AffineExpression(x)
convert(::Type{AffineExpression{T}}, x::Number) where T = AffineExpression(convert(T, x))
convert(::Type{AffineExpression{T}}, ae::AffineExpression) where {T<:Number} = AffineExpression{T}(ae.coeffs)
