"""
    eigenbox(A)

Returns an enclosure of all the eigenvalues of `A`. If `A` is symmetric, than the
output is a real interval, otherwise it is a complex interval.

### Algorithm

The algorithms used by the function are described in [[HLA13]](@ref).

### Notes

The enclosure is not rigorous, meaning that the real eigenvalue problems solved internally
utilize normal floating point computations.

### Examples

```jldoctest
julia> A = [0 -1 -1;2 -1.399.. -0.001 0;1 0.5 -1]
3×3 Matrix{Interval{Float64}}:
 [0, 0]  [-1, -1]                       [-1, -1]
 [2, 2]       [-1.39901, -0.000999999]    [0, 0]
 [1, 1]        [0.5, 0.5]               [-1, -1]

julia> eigenbox(A)
[-1.90679, 0.970154] + [-2.51903, 2.51903]im
```
"""
function eigenbox(A::Symmetric{Interval{T}, Matrix{Interval{T}}}) where {T}

    AΔ = Symmetric(radius.(A))
    Ac = Symmetric(mid.(A))

    ρ = eigmax(AΔ)
    λmax = eigmax(Ac)
    λmin = eigmin(Ac)
    return Interval(λmin - ρ, λmax + ρ)

end

function eigenbox(A::AbstractMatrix{Interval{T}}) where {T}

    λ = eigenbox(Symmetric(0.5*(A + A')))

    n = checksquare(A)
    μ = eigenbox(Symmetric([zeros(n, n) 0.5*(A - A');
                            0.5*(A' - A) zeros(n, n)]))

    return λ + μ*im
end

function eigenbox(M::Matrix{Complex{Interval{T}}}) where {T}
    A = real.(M)
    B = imag.(M)
    λ = eigenbox(Symmetric(0.5*[A+A' B'-B;
                                B-B' A+A']))

    μ = eigenbox(Symmetric(0.5*[B+B' A-A';
                                A'-A B+B']))

    return λ + μ*im
end


function eigenbox(M::Hermitian{Complex{Interval{T}}, Matrix{Complex{Interval{T}}}}) where T
    A = real(M)
    B = imag(M)
    return eigenbox(Symmetric([A B';B A]))

end
