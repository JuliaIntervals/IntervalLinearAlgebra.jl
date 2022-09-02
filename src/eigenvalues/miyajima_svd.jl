# In this file we implement some of the algorithms from
#
# Shinya Miyajima
#
# Verified bounds for all the singular values of matrix
# Japan J. Indust. Appl. Math.
# DOI 10.1007/s13160-014-0145-5 

abstract type AbstractIntervalSVDSolver end

struct M1 <: AbstractIntervalSVDSolver end

# The method is called svdbox to mirror eigenbox, even if it would be best to call it 
# svd interval

using LinearAlgebra

function _power_iteration_singularvalue(A, max_iter)
    n = size(A, 2)
    xp = rand(n)
    @inbounds for _ in 1:max_iter
    xp .= A'*(A*xp)
    xp ./= norm(xp)
    end
    return xp
end

# we use this function to bound the $2$ norm of a matrix,
# as remarked in Rump, Verified bounds for singular values,
# Equation 3.3
function _bound_perron_frobenius_singularvalue(M, max_iter=10)

    size(M) == (1, 1) && return M[1]
    xpf = IA.Interval.(_power_iteration_singularvalue(mid.(M), max_iter))
    Mxpf = M'* (M * xpf)
    ρ = zero(eltype(M))
    @inbounds for (i, xi) in enumerate(xpf)
        iszero(xi) && continue
        tmp = Mxpf[i] / xi
        ρ = max(ρ, tmp.hi)
    end
    return ρ
end


"""
    svdbox(A[, method = M1()])

Return an enclosure for all the singular values of `A`.

# Input 
-   `A` -- interval matrix
- `method` -- the method used to compute the enclosure
    Possible values are
        - M1 Algorithm from Theorem 7 in [[MIYA14]](@ref)

### Algorithm 

The algorithms used by the function are described in [[MIYA14]](@ref).

### Example
The matrix encloses a matrix that has singular values 
``3, √5, 2, 0``.

```julia
julia> A = [0.9..1.1 0 0 0 2; 0 0 3 0 0; 0 0 0 0 0; 0 2 0 0 0]
4×5 Matrix{Interval{Float64}}:
     [0.899999, 1.10001]  [0, 0]  [0, 0]  [0, 0]  [2, 2]
 [0, 0]                   [0, 0]  [3, 3]  [0, 0]  [0, 0]
 [0, 0]                   [0, 0]  [0, 0]  [0, 0]  [0, 0]
 [0, 0]                   [2, 2]  [0, 0]  [0, 0]  [0, 0]

 julia> IntervalLinearAlgebra.svdbox(A, IntervalLinearAlgebra.M1())
 4-element Vector{Interval{Float64}}:
   [2.98999, 3.01001]
   [2.22606, 2.24607]
   [1.98999, 2.01001]
  [-0.0100001, 0.0100001]
```

"""
function svdbox(A::AbstractMatrix{Interval{T}}, method = ::M1) where T
    mA = mid.(A)
    svdA  = svd(mA)
    U = Interval{T}.(svdA.U)
    Vt = Interval{T}.(svdA.Vt)
    Σ = diagm(Interval{T}.(svdA.S))
    V = Interval{T}.(svdA.V)

    E = U*Σ*Vt - A
    normE = sqrt(_bound_perron_frobenius_singularvalue(E))
    
    F = Vt*V-I
    normF = sqrt(_bound_perron_frobenius_singularvalue(F))

    G = U'*U-I
    normG = sqrt(_bound_perron_frobenius_singularvalue(G))

    @assert normF < 1 "It is not possible to verify the singular values with this precision"
    @assert normG < 1 "It is not possible to verify the singular values with this precision"

    svdbounds = [hull(σ*sqrt((1-normF)*(1-normG))-normE, 
                     σ*sqrt((1+normF)*(1+normG))+normE)
                for σ in diag(Σ)]

    return svdbounds
end