abstract type LinearSolver end


"""
Solves the linear system using Hansen-Bliek-Rohn method. This might give a tighter solution
than Gauss elimination in general, but this method only works for some matrices
(at least with H-matrices without preconditioning, I think it should work in general
for strongly regular matrices using preconditioning). See section 5.6.2 of [1] (page 49)
"""
struct HansenBliekRohn <: LinearSolver end

function (hbr::HansenBliekRohn)(A, b)
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

struct GaussElimination <: LinearSolver end

function (ge::GaussElimination)(A, b)
    n = length(b)
    A = MMatrix{n, n}(A)
    b = MVector{n}(b)
    # TODO: IMPLEMENT IT :D
end


## JACOBI
"""
Solves the linear system using Jacobi method. See section 5.7.4 of [1] (page 52)
"""
struct Jacobi <: LinearSolver
    max_iterations::Int
    atol::Float64
end

Jacobi() = Jacobi(20, 0.0)

function (jac::Jacobi)(x, A, b)

    n = length(b)
    for idx in 1:jac.max_iterations
        xold = copy(x)
        @inbounds @simd for i in 1:n
            x[i] = b[i]
            for j in 1:n
                (i == j) || (x[i] -= A[i, j] * xold[j])
            end
            x[i] = (x[i]/A[i, i]) ∩ xold[i]
        end
        all(isapprox.(x, xold; atol=jac.atol)) && break
    end
    nothing
end

## GAUSS SEIDEL
"""
Solves the linear system using Gauss-Seidel method. See section 5.7.4 of [1] (page 52)
"""
struct GaussSeidel <: LinearSolver
    max_iterations::Int
    atol::Float64
end

GaussSeidel() = GaussSeidel(20, 0.0)

function (gs::GaussSeidel)(x, A, b)
    n = length(b)

    @inbounds for _ in 1:gs.max_iterations
        xold = copy(x)
        @inbounds for i in 1:n
            x[i] = b[i]
            @inbounds for j in 1:n
                (i == j) || (x[i] -= A[i, j] * x[j])
            end
            x[i] = (x[i]/A[i, i]) .∩ xold[i]
        end
        all(isapprox.(x, xold; atol=gs.atol)) && break
    end
    nothing
end

## KRAWCZYK
"""
Solves the linear system using Krawczyk method. Note that Krawczyk does not work for
matrices whose middle matrix is diagonal. See section 5.7.3 of [1] (page 51)
"""
struct Krawczyk <: LinearSolver
    max_iterations::Int
    atol::Float64
end

Krawczyk() = Krawczyk(20, 0.0)

function (kra::Krawczyk)(x, A, b)
    C = inv(mid.(A))
    for i = 1:kra.max_iterations
        xnew  = (C*b  - C*(A*x) + x) .∩ x
        all(isapprox.(x, xnew; atol=kra.atol)) && return xnew
        x = xnew
    end
    return x
end


## wrapper

function solve(A, b, method)

    A, b = precondition(A, b)
    x = enclose(A, b)

    method(x, A, b)

    return x
end

function solve(A, b, method::HansenBliekRohn)
    A, b = precondition(A, b)
    return method(A, b)
end
