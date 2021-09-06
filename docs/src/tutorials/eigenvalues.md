# Eigenvalue computations

## Eigenvalues of interval matrices

Given a (real or complex) interval matrix ``A\in\mathbb{IC}^{n\times n}``, we define the eigenvalue set 

```math
\mathbf{\Lambda}=\{\lambda\in\mathbb{C}: \lambda\text{ is an eigenvalue of }A\text{ for some }A\in\mathbf{A}\}.
```

While characterizing the solution set ``\mathbf{\Lambda}`` (or even its hull) is computationally challenging, the package offers the function [`eigenbox`](@ref) which contains an interval box containing ``\mathbf{\Lambda}``. 

!!! note
    At the moment, `eigenbox` is not rigorous, that is the computations for the non-interval eigenvalue problem solved internally are carried out using normal non-verified floating point computations.

To demonstrate the functionality, let us consider the following interval matrix

```@example eigs
using IntervalLinearAlgebra

A = [-3.. -2 4..5 4..6 -1..1.5;
    -4.. -3 -4.. -3 -4.. -3 1..2;
    -5.. -4 2..3 -5.. -4 -1..0;
    -1..0.1 0..1 1..2 -4..2.5]
```

Now we can bound the eigenvalue set
```@example eigs
ebox = eigenbox(A)
```

To get a qualitative evaluation of the enclosure, we can simulate the solution set of ``\mathbf{A}`` using Montecarlo, as it is done in the following example

```@example eigs
using Plots
N = 1000

evalues = zeros(ComplexF64, 4, N)

for i in 1:N
    evalues[:, i] = eigvals(rand.(A))
end

rpart = real.(evalues)
ipart = imag.(evalues)

plot(IntervalBox(real(ebox), imag(ebox)); ratio=1, label="enclosure")
scatter!(rpart[1, :], ipart[1, :]; label="λ₁")
scatter!(rpart[2, :], ipart[2, :]; label="λ₂")
scatter!(rpart[3, :], ipart[3, :]; label="λ₃")
scatter!(rpart[4, :], ipart[4, :]; label="λ₄")
xlabel!("real")
ylabel!("imag")
savefig("eigs.png") # hide
```

![](eigs.png)

Internally, the generical interval eigenvalue proper is reduced to a symmetric interval eigenvalue problem, as described in [[HLA13]](@ref). It is good to highlight that The symmetric interval eigenvalue problem can be solved in two ways

- Rohn method -- (default one) computes an enclosure of the eigenvalues set for the symmetric interval matrix. This is fast but the enclosure can be strictly larger than the hull
- Hertz method -- computes the exact hull of the eigenvalues for the symmetric interval matrix. Generally, these leads to tigher bounds, but it has exponential complexity, so it will be unfeasible for big matrices.

The function `eigenbox` takes a second optional parameter (RohnMethod by default) to specify what algorithm to use for the symmetric interval eigenvalue problem. The following example bounds the eigenvalues of the previous matrix using HertzMethod, as can be noticed by the figure below, the Hertz method gives a tighter bound on the eigenvalues set.

```@example eigs
eboxhertz = eigenbox(A, HertzMethod)
```

```@example eigs
plot(IntervalBox(real(ebox), imag(ebox)); ratio=1, label="enclosure")
plot!(IntervalBox(real(eboxhertz), imag(eboxhertz)); label="Hertz enclosure", color="#00FF00") # hide
scatter!(rpart[1, :], ipart[1, :]; label="λ₁") # hide
scatter!(rpart[2, :], ipart[2, :]; label="λ₂") # hide
scatter!(rpart[3, :], ipart[3, :]; label="λ₃") # hide
scatter!(rpart[4, :], ipart[4, :]; label="λ₄") # hide
xlabel!("real")
ylabel!("imag")
savefig("eigs2.png") # hide
```

![](eigs2.png)

## Verified floating point computations of eigenvalues

In the previous section we considered the problem of finding the solution set (or an enclosure of it) of an interval matrix. In this section, we consider the problem of computing eigenvalues and eigenvectors of a floating point matrix *rigorously*, that is we want to find an enclosure of the true eigenvalues and eigenvectors of the matrix. In `IntervalLinearAlgebra.jl` this is achieved using the [`verify_eigen`](@ref) function, as the following example demonstrates.

```@example eigs
A = [1 2; 3 4]
evals, evecs, cert = verify_eigen(A)
evals
```

```@example eigs
evecs
```

```@example eigs
cert
```

If called with only one input `verify_eigen` will first compute an approximate solution for the eigenvalues and eigenvectors of ``A`` and use that to find a rigorous bounding on the true eigenvalues and eigenvectors of the matrix. It is also possible to give the function the scalar parameters ``\lambda`` and ``\vec{v}``, in which case it will compute a rigorous bound only for the specified eigenvalue ``\lambda`` and eigenvector ``\vec{v}``. The last output of the function is a vector of boolean certificates, if the ``i``th element is set to true, then the enclosure of the ``i``th eigenvalue and eigenvector is rigorous, that is the algorithm could prove that that enclosure contains the true eigenvalue and eigenvector of ``A``. If the certificate is false, then the algorithm could not prove the validity of the enclosure.

The function also accepts interval inputs. This is handy if the input matrix elements cannot be represented exactly as floating point numbers. Note however that this is meant only for interval matrices with very small intervals. If you have larger intervals, you should use the function of the previous section.

To test the function, let us consider the following example. First we generate random eigenvalues and eigenvectors

```@example eigs
using Random; # hide
Random.seed!(42) # hide
ev = sort(randn(10))
D = Diagonal(ev)
P = randn(10, 10)
Pinv, _ = epsilon_inflation(P, Diagonal(ones(10)))
A = interval.(P) * D * Pinv
```

Now we obtained an interval matrix ``\mathbf{A}`` so that `ev` and `P` are eigenvalues and eigenvectors of some ``A\in\mathbf{A}``. Note that ``P^{-1}`` had to be computed rigorously using [`epsilon_inflation`](@ref). Now we can compute its eigenvalues and eigenvectors and verify that the enclosures contain the true values.

```@example eigs
evals, evecs, cert = verify_eigen(A)
evals
```

```@example eigs
evecs
```

```@example eigs
cert
```

```@example eigs
ev .∈ evals
```

Note also that despite the original eigenvalues and eigenvectors were real, the returned enclosures are complex. This is because any infinitesimally small perturbation in the elements of ``A`` may cause the eigenvalues to move away from the real line. For this reason, unless the matrix has some special structure that guarantees the eigenvalues are real (e.g. symmetric matrices), a valid enclosure should always be complex.

Finally, the concept of enclosure of eigenvector may feel confusing, since eigenvectors are unique up to scale.
This scale ambiguity is resolved by starting with the approximate eigenvector computed by normal linear algebra routines and fixing the element with the highest magnitude. 

```@example eigs
complex_diam(x) = max(diam(real(x)), diam(imag(x)))

complex_diam.(evecs)
```

As can be seen, for each eigenvector there's an interval with zero width, since to resolve scale ambiguity one non-zero element can be freely chosen (assuming eigenvalues have algebraic multiplicity ``1``). After that, the eigenvector is fixed and it makes sense to talk about enclosures of the other elements.


