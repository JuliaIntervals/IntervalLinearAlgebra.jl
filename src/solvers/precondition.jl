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
    InverseMidpoint <: Precondition

Preconditioner that preconditions the linear system ``Ax=b`` with ``inv(Ac)``, where `Ac` is the midpoint matrix of ``A``.

### Example
```julia
imp = InverseMidpoint() # instantiate preconditioner
Aprec, bprec = imp(A, b) # apply preconditioner
```
"""
struct InverseMidpoint <: Precondition end

function (icp::InverseMidpoint)(A, b)
    R = inv(mid.(A))
    return R*A, R*b
end

"""
    InverseDiagonalMidpoint <: Precondition

Preconditioner that preconditions the linear system ``Ax=b`` with ``inv(Diagonal(Ac))``, where `Ac` is the midpoint matrix of ``A``.

### Example
```julia
idmp = InverseDiagonalMidpoint() # instantiate preconditioner
Aprec, bprec = idmp(A, b) # apply preconditioner
```
"""
struct InverseDiagonalMidpoint <: Precondition end

function (idp::InverseDiagonalMidpoint)(A, b)
    R = inv(Diagonal(mid.(A)))
    return R*A, R*b
end