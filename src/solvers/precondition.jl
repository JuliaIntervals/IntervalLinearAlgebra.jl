abstract type Precondition end

"""
    NoPrecondition <: Precondition

Trivial preconditioner which does nothing.

### Example
```julia
np = NoPrecondition() # instantiate preconditioner
Aprec, bprec = np(A, b) # apply preconditioner
```
"""
struct NoPrecondition <: Precondition end

(np::NoPrecondition)(A, b) = A, b

"""
    InverseMidpointPrecondition <: Precondition

Preconditioner that preconditions the linear system ``Ax=b`` with ``inv(Ac)``, where `Ac` is the midpoint matrix of ``A``.

### Example
```julia
imp = InverseMidpointPrecondition() # instantiate preconditioner
Aprec, bprec = imp(A, b) # apply preconditioner
```
"""
struct InverseMidpointPrecondition <: Precondition end

function (icp::InverseMidpointPrecondition)(A, b)
    R = inv(mid.(A))
    return R*A, R*b
end

"""
    InverseDiagonalMidpointPrecondition <: Precondition

Preconditioner that preconditions the linear system ``Ax=b`` with ``inv(Diagonal(Ac))``, where `Ac` is the midpoint matrix of ``A``.

### Example
```julia
idmp = InverseDiagonalMidpointPrecondition() # instantiate preconditioner
Aprec, bprec = idmp(A, b) # apply preconditioner
```
"""
struct InverseDiagonalMidpointPrecondition <: Precondition end

function (idp::InverseDiagonalMidpointPrecondition)(A, b)
    R = inv(Diagonal(mid.(A)))
    return R*A, R*b
end