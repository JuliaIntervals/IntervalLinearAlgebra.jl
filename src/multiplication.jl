const config = Dict(:multiplication => :fast)

struct MultiplicationType{T} end

getconfiguration() = config

"""
    setmultiplication(multype)

Sets the algorithm used to perform matrix multiplication with interval matrices.

### Input

- `multype` -- symbol describing the algorithm used
   - `:slow` -- uses traditional matrix multiplication algorithm.
   - `:fast` -- computes an enclosure of the matrix product using the midpoint-radius
                notation of the matrix [`RUM10`](@ref).

### Notes

- By default, `:fast` is used.
- Using `fast` is generally significantly faster, but it may return larger intervals
  (50% wider at most).
"""
function setmultiplication(multype)
    type = MultiplicationType{multype}()
    @eval *(A::AbstractMatrix{Interval{T}} where T, B::AbstractMatrix{Interval{T}} where T) =
        *($type, A, B)

    config[:multiplication] = multype
end

function *(::MultiplicationType{:slow}, A, B)
    TS = promote_type(eltype(A), eltype(B))
    return mul!(similar(B, TS, (size(A,1), size(B,2))), A, B)
end

function *(::MultiplicationType{:fast},
           A::AbstractMatrix{Interval{T}},
           B::AbstractMatrix{Interval{T}}) where {T<:Real}

    Ainf = inf.(A)
    Asup = sup.(A)
    Binf = inf.(B)
    Bsup = sup.(B)

    mA, mB, R, Csup = setrounding(T, RoundUp) do
        mA = Ainf + 0.5 * (Asup - Ainf)
        mB = Binf + 0.5 * (Bsup - Binf)

        rA = mA - Ainf
        rB = mB - Binf

        R = abs.(mA) * rB + rA * (abs.(mB) + rB)
        Csup = mA * mB + R

        return mA, mB, R, Csup
    end

    Cinf = setrounding(T, RoundDown) do
        mA * mB - R
    end

    return Interval.(Cinf, Csup)
end
