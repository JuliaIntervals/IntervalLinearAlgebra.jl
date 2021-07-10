"""
    AbstractLinearSolver

Abstract type for solvers of interval linear systems.
"""
abstract type AbstractLinearSolver end

"""
    AbstractDirectSolver <: AbstractLinearSolver

Abstract type for direct solvers of interval linear systems, such as Gaussian elimination
and Hansen-Bliek-Rohn.
"""
abstract type AbstractDirectSolver <: AbstractLinearSolver end

"""
    AbstractIterativeSolver <: AbstractLinearSolver

Abstract type for iterative solvers of interval linear systems, such as Jacobi or
Gauss-Seidel.
"""
abstract type AbstractIterativeSolver <: AbstractLinearSolver end

"""
    HansenBliekRohn <: AbstractDirectSolver

Type for the `HansenBliekRohn` solver of the square interval linear system ``Ax=b``.
For more details see section 5.6.2 of [[HOR19]](@ref)

### Notes

- Hansen-Bliek-Rohn works with H-matrices without precondition and with strongly regular
  matrices using [`InverseMidpoint`](@ref) precondition
- If the midpoint of ``A`` is a diagonal matrix, then the algorithm returns the exact hull.
- An object of type Hansen-Bliek-Rohn is a callable function with method

        (hbr::HansenBliekRohn)(A::AbstractMatrix{T},
                               b::AbstractVector{T}) where {T<:Interval}

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-1, 1]

julia> hbr = HansenBliekRohn()
HansenBliekRohn linear solver

julia> hbr(A, b)
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]
 [-1.33334, 1.33334]
```
"""
struct HansenBliekRohn <: AbstractDirectSolver end

function (hbr::HansenBliekRohn)(A::AbstractMatrix{T},
                                b::AbstractVector{T}) where {T<:Interval}
    n = length(b)
    compA = comparison_matrix(A)
    compA_inv = inv(compA)
    u = compA_inv*mag.(b)
    d = diag(compA_inv)
    α = diag(compA) .- 1 ./d
    α = Interval.(-α, α) #TODO: probably directed rounded is needed here, need to check
    β = @. u/d - mag(b)
    β = Interval.(-β, β)
    x = (b .+ β)./(diag(A) .+ α)

end

"""
    GaussElimination <: AbstractDirectSolver

Type for the Gaussian elimination solver of the square interval linear system ``Ax=b``.
For more details see section 5.6.1 of [[HOR19]](@ref)

### Notes

- An object of type `GaussianElimination` is a callable function with method

        (ge::GaussianElimination)(A::AbstractMatrix{T},
                                  b::AbstractVector{T}) where {T<:Interval}

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-1, 1]

julia> ge = GaussianElimination()
GaussianElimination linear solver

julia> ge(A, b)
2-element Vector{Interval{Float64}}:
 [-1.66667, 1.66667]
 [-1.33334, 1.33334]
```
"""
struct GaussianElimination <: AbstractDirectSolver end

function (ge::GaussianElimination)(A::AbstractMatrix{T},
                                   b::AbstractVector{T}) where {T<:Interval}
    n = length(b)
    Abrref = rref([A b])

    # backsubstitution
    x = similar(b)
    x[end] = Abrref[n, n+1]/Abrref[n, n]
    @inbounds for i = n-1:-1:1
        x[i] = (Abrref[i, n+1] - sum(Abrref[i, j]*x[j] for j in i+1:n))/Abrref[i, i]
    end
    return x
end


## JACOBI
"""
    Jacobi <: AbstractIterativeSolver

Type for the Jacobi solver of the interval linear system ``Ax=b``.
For details see Section 5.7.4 of [[HOR19]](@ref)

### Fields

- `max_iterations` -- maximum number of iterations (default 20)
- `atol` -- absolute tolerance (default 0), if at some point ``|xₖ - xₖ₊₁| < atol``
            (elementwise), then stop and return ``xₖ₊₁``.
            If `atol=0`, then `min(diam(A))*1e-5` is used.

### Notes

- An object of type `Jacobi` is a function with method

        (jac::Jacobi)(A::AbstractMatrix{T},
                      b::AbstractVector{T},
                      [x]::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

    #### Input
    - `A`   -- N×N interval matrix
    - `b`   -- interval vector of length N
    - `x`   -- (optional) initial enclosure for the solution of ``Ax = b``. If not given,
                it is automatically computed using [`enclose`](@ref enclose)

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-1, 1]

julia> jac = Jacobi()
Jacobi linear solver
max_iterations = 20
atol = 0.0

julia> jac(A, b)
2-element Vector{Interval{Float64}}:
 [-1.66668, 1.66668]
 [-1.33335, 1.33335]
```
"""
struct Jacobi <: AbstractIterativeSolver
    max_iterations::Int
    atol::Float64
end

Jacobi() = Jacobi(20, 0.0)

function (jac::Jacobi)(A::AbstractMatrix{T},
                       b::AbstractVector{T},
                       x::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

    n = length(b)
    atol = iszero(jac.atol) ? minimum(diam.(A))*1e-5 : jac.atol

    for _ in 1:jac.max_iterations
        xold = copy(x)
        @inbounds @simd for i in 1:n
            x[i] = b[i]
            for j in 1:n
                (i == j) || (x[i] -= A[i, j] * xold[j])
            end
            x[i] = (x[i]/A[i, i]) ∩ xold[i]
        end
        all(interval_isapprox.(x, xold; atol=atol)) && break
    end
    return x
end

## GAUSS SEIDEL
"""
    GaussSeidel <: AbstractIterativeSolver

Type for the Gauss-Seidel solver of the interval linear system ``Ax=b``.
For details see Section 5.7.4 of [[HOR19]](@ref)

### Fields

- `max_iterations` -- maximum number of iterations (default 20)

- `atol` -- absolute tolerance (default 0), if at some point ``|xₖ - xₖ₊₁| < atol``
            (elementwise), then stop and return ``xₖ₊₁``.
            If `atol=0`, then `min(diam(A))*1e-5` is used.

### Notes

- An object of type `GaussSeidel` is a function with method

        (gs::GaussSeidel)(A::AbstractMatrix{T},
                          b::AbstractVector{T},
                          [x]::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

    #### Input
    - `A`   -- N×N interval matrix
    - `b`   -- interval vector of length N
    - `x`   -- (optional) initial enclosure for the solution of ``Ax = b``. If not given,
                it is automatically computed using [`enclose`](@ref enclose)

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-1, 1]

julia> gs = GaussSeidel()
GaussSeidel linear solver
max_iterations = 20
atol = 0.0

julia> gs(A, b)
2-element Vector{Interval{Float64}}:
 [-1.66668, 1.66668]
 [-1.33334, 1.33334]
```
"""
struct GaussSeidel <: AbstractIterativeSolver
    max_iterations::Int
    atol::Float64
end

GaussSeidel() = GaussSeidel(20, 0.0)

function (gs::GaussSeidel)(A::AbstractMatrix{T},
                           b::AbstractVector{T},
                           x::AbstractVector{T}=enclose(A, b)) where {T<:Interval}
    n = length(b)

    atol = iszero(gs.atol) ? minimum(diam.(A))*1e-5 : gs.atol
    @inbounds for _ in 1:gs.max_iterations
        xold = copy(x)
        @inbounds for i in 1:n
            x[i] = b[i]
            @inbounds for j in 1:n
                (i == j) || (x[i] -= A[i, j] * x[j])
            end
            x[i] = (x[i]/A[i, i]) .∩ xold[i]
        end
        all(interval_isapprox.(x, xold; atol=atol)) && break
    end
    return x
end

## KRAWCZYK
"""
    LinearKrawczyk <: AbstractIterativeSolver

Type for the Krawczyk solver of the interval linear system ``Ax=b``.
For details see Section 5.7.3 of [[HOR19]](@ref)

### Fields

- `max_iterations` -- maximum number of iterations (default 20)
- `atol` -- absolute tolerance (default 0), if at some point ``|xₖ - xₖ₊₁| < atol``
            (elementwise), then stop and return ``xₖ₊₁``.
            If `atol=0`, then `min(diam(A))*1e-5` is used.

### Notes

- An object of type `LinearKrawczyk` is a function with method

        (kra::LinearKrawczyk)(A::AbstractMatrix{T},
                              b::AbstractVector{T},
                              [x]::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

    #### Input
    - `A`   -- N×N interval matrix
    - `b`   -- interval vector of length N
    - `x`   -- (optional) initial enclosure for the solution of ``Ax = b``. If not given,
                it is automatically computed using [`enclose`](@ref enclose)

### Examples

```jldoctest
julia> A = [2..4 -1..1;-1..1 2..4]
2×2 Matrix{Interval{Float64}}:
  [2, 4]  [-1, 1]
 [-1, 1]   [2, 4]

julia> b = [-2..2, -1..1]
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-1, 1]

julia> kra = LinearKrawczyk()
LinearKrawczyk linear solver
max_iterations = 20
atol = 0.0


julia> kra(A, b)
2-element Vector{Interval{Float64}}:
 [-2, 2]
 [-2, 2]
```
"""
struct LinearKrawczyk <: AbstractIterativeSolver
    max_iterations::Int
    atol::Float64
end

LinearKrawczyk() = LinearKrawczyk(20, 0.0)

function (kra::LinearKrawczyk)(A::AbstractMatrix{T},
                               b::AbstractVector{T},
                               x::AbstractVector{T}=enclose(A, b)) where {T<:Interval}

    atol = iszero(kra.atol) ? minimum(diam.(A))*1e-5 : kra.atol

    C = inv(mid.(A))
    for i = 1:kra.max_iterations
        xnew  = (C*b  - C*(A*x) + x) .∩ x
        all(interval_isapprox.(x, xnew; atol=atol)) && return xnew
        x = xnew
    end
    return x
end

# custom printing for solvers
function Base.string(s::AbstractLinearSolver)

    str="""$(typeof(s)) linear solver
    """

    fields = fieldnames(typeof(s))
    for field in fields
        str *= """$field = $(getfield(s, field))
        """
    end
    return str
end

Base.show(io::IO, s::AbstractLinearSolver) = print(io, string(s))
