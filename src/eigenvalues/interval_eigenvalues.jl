abstract type AbstractIntervalEigenSolver end

struct Hertz <: AbstractIntervalEigenSolver end
struct Rohn <: AbstractIntervalEigenSolver end

"""
    eigenbox(A[, method=Rohn()])

Returns an enclosure of all the eigenvalues of `A`. If `A` is symmetric, then the
output is a real interval, otherwise it is a complex interval.

### Input

- `A` -- square interval matrix
- `method` -- method used to solve the symmetric interval eigenvalue problem (bounding
    eigenvalues of general matrices is also reduced to the symmetric case).
    Possible values are

      - `Rohn` -- (default) fast method to compute an enclosure of the eigenvalues of
            a symmetric interval matrix
      - `Hertz` -- finds the exact hull of the eigenvalues of a symmetric interval
            matrix, but has exponential complexity.

### Algorithm

The algorithms used by the function are described in [[HLA13]](@ref).

### Notes

The enclosure is not rigorous, meaning that the real eigenvalue problems solved internally
utilize normal floating point computations.

### Examples

```jldoctest
julia> A = [0 -1 -1; 2 -1.399.. -0.001 0; 1 0.5 -1]
3×3 Matrix{Interval{Float64}}:
 [0.0, 0.0]_com_NG  [-1.0, -1.0]_com_NG   [-1.0, -1.0]_com_NG
 [2.0, 2.0]_com_NG  [-1.399, -0.001]_com   [0.0, 0.0]_com_NG
 [1.0, 1.0]_com_NG   [0.5, 0.5]_com_NG    [-1.0, -1.0]_com_NG

julia> eigenbox(A)
[-1.90678, 0.970154]_com_NG + im*[-2.51903, 2.51903]_com_NG

julia> eigenbox(A, Hertz())
[-1.64731, 0.520455]_com_NG + im*[-2.11119, 2.11119]_com_NG
```
"""
function eigenbox(A::Symmetric{Interval{T}, Matrix{Interval{T}}}, ::Rohn) where {T}

    AΔ = Symmetric(IntervalArithmetic.radius.(A))
    Ac = Symmetric(mid.(A))

    ρ = eigmax(AΔ)
    λmax = eigmax(Ac)
    λmin = eigmin(Ac)
    return IA.interval(λmin - ρ, λmax + ρ)

end


function eigenbox(A::Symmetric{Interval{T}, Matrix{Interval{T}}}, ::Hertz) where {T}

    n = checksquare(A)
    Amax = Matrix{T}(undef, n, n)
    Amin = Matrix{T}(undef, n, n)

    λmin = Inf
    λmax = -Inf
    @inbounds for z in Orthants(n)
        first(z) < 0 && continue
        for j in 1:n
            for i in 1:j
                if z[i] == z[j]
                    Amax[i, j] = sup(A[i, j])
                    Amin[i, j] = inf(A[i, j])
                else
                    Amax[i, j] = inf(A[i, j])
                    Amin[i, j] = sup(A[i, j])
                end
            end
        end

        candmax = eigmax(Symmetric(Amax))
        candmin = eigmin(Symmetric(Amin))
        λmin = min(λmin, candmin)
        λmax = max(λmax, candmax)
    end
    return IA.interval(λmin, λmax)
end

function eigenbox(A::AbstractMatrix{Interval{T}},
                  method::AbstractIntervalEigenSolver) where {T}

    λ = eigenbox(Symmetric(0.5*(A + A')), method)

    n = size(A, 1)
    μ = eigenbox(Symmetric([zeros(n, n) 0.5*(A - A');
                            0.5*(A' - A) zeros(n, n)]), method)

    return λ + μ*im
end

function eigenbox(M::AbstractMatrix{Complex{Interval{T}}},
                  method::AbstractIntervalEigenSolver) where {T}
    A = real.(M)
    B = imag.(M)
    λ = eigenbox(Symmetric(0.5*[A+A' B'-B;
                                B-B' A+A']), method)

    μ = eigenbox(Symmetric(0.5*[B+B' A-A';
                                A'-A B+B']), method)

    return λ + μ*im
end


function eigenbox(M::Hermitian{Complex{Interval{T}}, Matrix{Complex{Interval{T}}}},
                  method::AbstractIntervalEigenSolver) where {T}
    A = real(M)
    B = imag(M)
    return eigenbox(Symmetric([A B';B A]), method)
end

# default
eigenbox(A) = eigenbox(A, Rohn())
