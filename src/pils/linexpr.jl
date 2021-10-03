import Base: +, -, *, ==, show

struct AffineExpression{T}
    coeffs::Vector{T}
end

## VARIABLES CONSTRUCTION

const _vars_dict = Dict(:vars => Symbol[])

macro linvars(x...)
    _vars_dict[:vars] = Symbol[]
    if length(x) == 1
        s = x[1].args[1]
        start = x[1].args[2].args[2]
        stop = x[1].args[2].args[3]
        _vars_dict[:vars] = [Symbol(s, i) for i in start:stop]
    else
        _vars_dict[:vars] = collect(x)
    end

    n = length(_vars_dict[:vars])
    ex = quote end

    for (i, s) in enumerate(_vars_dict[:vars])
        c = zeros(Int, n + 1)
        c[i] = 1
        push!(ex.args, :($(esc(s)) = AffineExpression($c)))
    end
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
# ARITHMETIC OPERATIONS #
#########################

for op in (:+, :-)
    @eval function $op(ae1::AffineExpression, ae2::AffineExpression)
        return AffineExpression($op(ae1.coeffs, ae2.coeffs))
    end

    @eval function $op(ae::AffineExpression, n::Number)
        cnew = copy(ae.coeffs)
        cnew[end] = $op(cnew[end], n)
        return AffineExpression(cnew)
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
## CONVERTIONS
