"""
    verify_eigen(A[, λ, X0]; w=0.1, ϵ=1e-16, maxiter=10)

Finds a rigorous bound for the eigenvalues and eigenvectors of `A`. Eigenvalues are treated
as simple.

### Input

- `A`       -- matrix
- `λ`       -- (optional) approximate value for an eigenvalue of `A`
- `X0`      -- (optional) eigenvector associated to `λ`
- `w`       -- relative inflation parameter
- `ϵ`       -- absolute inflation parameter
- `maxiter` -- maximum number of iterations

### Output

- Interval bounds on eigenvalues and eigenvectors.
- A boolean certificate (or a vector of booleans if all eigenvalues are computed) `cert`.
  If `cert[i]==true`, then the bounds for the ith eigenvalue and eigenvectore are rigorous,
  otherwise not.

### Algorithm

The algorithm for this function is described in [[RUM01]](@ref).

### Example

```julia
julia> A = Symmetric([1 2;2 3])
2×2 Symmetric{Int64, Matrix{Int64}}:
 1  2
 2  3

julia> evals, evecs, cert = verify_eigen(A);

julia> evals
2-element Vector{Interval{Float64}}:
 [-0.236068, -0.236067]
  [4.23606, 4.23607]

julia> evecs
2×2 Matrix{Interval{Float64}}:
 [-0.850651, -0.85065]  [0.525731, 0.525732]
  [0.525731, 0.525732]  [0.85065, 0.850651]

julia> cert
2-element Vector{Bool}:
 1
 1
```
"""
function verify_eigen(A; kwargs...)
    evals, evecs = eigen(mid.(A))

    T = interval_eigtype(A, evals[1])
    evalues = similar(evals, T)
    evectors = similar(evecs, T)

    cert = Vector{Bool}(undef, length(evals))
    @inbounds for (i, λ₀) in enumerate(evals)
        λ, v, flag = verify_eigen(A, λ₀, view(evecs, :,i); kwargs...)
        evalues[i] = λ
        evectors[:, i] .= v
        cert[i] = flag

    end
    return evalues, evectors, cert
end

function verify_eigen(A, λ, X0; kwargs...)
    ρ, X, cert = _verify_eigen(A, λ, X0; kwargs...)
    return (real(λ) ± ρ) + (imag(λ) ± ρ) * im, X0 + X, cert
end

function verify_eigen(A::Symmetric, λ, X0; kwargs...)
    ρ, X, cert = _verify_eigen(A, λ, X0; kwargs...)
    return λ ± ρ, X0 + X, cert
end

function _verify_eigen(A, λ::Number, X0::AbstractVector;
                      w=0.1, ϵ=floatmin(), maxiter=10)

    _, v = findmax(abs.(X0))

    R = mid.(A) - λ * I
    R[:, v] .= -X0
    R = inv(R)
    C = IA.Interval.(A) - λ * I
    Z = -R * (C * X0)
    C[:, v] .= -X0
    C = I - R * C
    Zinfl = w * IA.Interval.(-mag.(Z), mag.(Z)) .+ IA.Interval(-ϵ, ϵ)

    X = Z
    cert = false
    @inbounds for _ in 1:maxiter
        if eltype(X) <: Complex # TODO: use dispatch
            Y = (real.(X) + Zinfl) + (imag.(X) + Zinfl) * im
        else
            Y = X + Zinfl
        end

        Ytmp = Y * Y[v]
        Ytmp[v] = 0

        X = Z + C * Y + R * Ytmp
        cert = all(X .⊂ Y)
        cert && break
    end

    ρ = mag(X[v])
    X[v] = 0

    return ρ, X, cert
end


"""
    bound_perron_frobenius_eigenvalue(A, max_iter=10)

Finds an upper bound for the Perron-Frobenius eigenvalue of the **non-negative** matrix `A`.

### Input

- `A` -- square real non-negative matrix
- `max_iter` -- maximum number of iterations of the power method used internally to compute
    an initial approximation of the Perron-Frobenius eigenvector

### Example

```julia-repl
julia> A = [1 2;3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> bound_perron_frobenius_eigenvalue(A)
5.372281323275249
```
"""
function bound_perron_frobenius_eigenvalue(A::AbstractMatrix{T}, max_iter=10) where {T<:Real}
    any(A .< 0) && throw(ArgumentError("Matrix contains negative entries"))
    return _bound_perron_frobenius_eigenvalue(A, max_iter)
end

function _bound_perron_frobenius_eigenvalue(M, max_iter=10)

    size(M, 1) == 1 && return M[1]
    xpf = IA.Interval.(_power_iteration(M, max_iter))
    Mxpf = M * xpf
    ρ = zero(eltype(M))
    @inbounds for (i, xi) in enumerate(xpf)
        iszero(xi) && continue
        tmp = Mxpf[i] / xi
        ρ = max(ρ, tmp.hi)
    end
    return ρ
end

function _power_iteration(A, max_iter)
    n = size(A,1)
    xp = rand(n)
    @inbounds for _ in 1:max_iter
    xp .= A*xp
    xp ./= norm(xp)
    end
    return xp
end


interval_eigtype(::Symmetric, ::T) where {T<:Real} = Interval{T}
interval_eigtype(::AbstractMatrix, ::T) where {T<:Real} = Complex{Interval{T}}
interval_eigtype(::AbstractMatrix, ::Complex{T}) where {T<:Real} = Complex{Interval{T}}
