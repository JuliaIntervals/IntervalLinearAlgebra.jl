"""
    verify_eigenvalues(A, λ, X0; k=0.1, ϵ=1e-16, maxiter=10)

Finds a rigorous bound for a given eigenvalue and the associated eigenvector(s) starting
from the approximate solutions `λ`, `X0`.

### Input

- `A`  -- matrix whose eigenvalue we want to compute
- `λ`  -- approximate value for an eigenvalue of `A`
- `X0` -- eigenvectors associated to `λ`, one for each column of X0

### Output

- An interval bound for `λ`
- An interval bound for `X0`
- A boolean certificate `cert`. If `cert==true`, then the bounds are rigorous, otherwise
  not.
"""
function verify_eigenvalues(A::Symmetric{T, Matrix{T}},
                            λ::Number, X0;
                            k=0.1, ϵ=floatmin(), maxiter=10) where {T}
    n = size(X0, 1)
    k = size(X0, 2)

    u = 1:n-k
    v = n-k+1:n
    R = mid.(A) - λ * I
    R[:, v] = -X0
    R = inv(R)
    C = IA.Interval.(A) - λ * I
    Z = -R * (C * X0)
    C[:, v] = -X0
    C = I - R * C

    Zinfl = IA.Interval.(-0.1*mag.(Z) .- ϵ, 0.1*mag.(Z) .+ ϵ)
    X = Z
    cert = false
    for _ in 1:maxiter
        Y = X + Zinfl
        Ytmp = copy(Y)
        Ytmp[v, :] .= 0
        X = Z + C * Y + R * (Ytmp * Y[v, :])
        cert = all(X .⊂ Y)
        cert && break
    end

    M = mag.(view(X, v, :))
    @show length(M)
    ρ = _bound_perron_frobenius_eigenvalue(M)

    X[v, :] .= 0
    return λ ± ρ, X0 + X, cert
end

"""
    verify_eigenvalues(A; k=0.1, ϵ=1e-16, maxiter=10)

Finds a rigorous bound for the eigenvalues and eigenvectors of `A`. Eigenvalues are treated
as simple.

### Input

- `A`  -- matrix whose eigenvalue we want to compute

### Output

- An object of type `Eigen` containing interval eigenvalues and eigenvectors.
- A boolean vector of certificates `cert`. If `cert[i]==true`, then the bounds for the ith
eigenvalue and eigenvectore are rigorous, otherwise not.

### Algorithm

The algorithm for this function is described in [[RUM99b]](@ref).

### Example

```julia
julia> A = [1 2;3 4]
2×2 Matrix{Int64}:
    1  2
    3  4

julia> ev, cert = verify_eigenvalues(A);

julia> ev
Eigen{Interval{Float64}, Interval{Float64}, Matrix{Interval{Float64}}, Vector{Interval{Float64}}}
values:
2-element Vector{Interval{Float64}}:
    [-0.372282, -0.372281]
    [5.37228, 5.37229]
vectors:
2×2 Matrix{Interval{Float64}}:
    [-0.824565, -0.824564]  [-0.415974, -0.415973]
    [0.565767, 0.565768]   [-0.909377, -0.909376]

julia> cert
2-element Vector{Bool}:
    1
    1
```
"""
function verify_eigenvalues(A::Symmetric{T, Matrix{T}}; kwargs...) where {T<:Real}
    tmp = eigen(mid.(A))
    ev = IA.Interval.(similar(tmp.values))
    vectors = IA.Interval.(similar(tmp.vectors))
    cert = similar(tmp.values, Bool)
    for i in 1:length(ev)
        λ, v, flag = verify_eigenvalues(A, tmp.values[i], tmp.vectors[:,i]; kwargs...)
        ev[i] = λ
        vectors[:, i] = v
        cert[i] = flag

    end
    return Eigen(ev, vectors), cert
end
