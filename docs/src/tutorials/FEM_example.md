# A simple FEM example

```@contents
Pages = ["FEM_example.md"]
```

In this section, a problem based on Example 4.1 from [https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/files/7271616/skalna2006.pdf] is considered. The matrix of the system is obtained using the Finite Element Method.

The stiffness matrix of a truss element in the local coordinate system is given by
```math
K_L = \frac{E A}{L}
\left(
  \begin{matrix}
  1 & 0 & -1 & 0 \\
  0 & 0 &  0 & 0 \\
 -1 & 0 &  1 & 0 \\
  0 & 0 &  0 & 0
\right),
```
and the change-of-basis matrix is given by
```math
_G(Q)_L = Q =
\left(
  \begin{matrix}
  cos(alpha) & -sin(alpha) & 0 & 0 \\
  sin(alpha) & cos(alpha) &  0 & 0 \\
  0 & 0 &  cos(alpha) & sin(alpha) \\
  0 & 0 &  sin(alpha) & cos(alpha)
\right),
```

The system of equations for each element is written in local coordinates as
```math
K_L d_L = f_L
```
and using the change-of-basis we obtain
```math
K_G d_G = f_G \qquad K_G = Q^T K_L Q
```

For element 1->3 the global stiffness matrix is
```math
K_1 = s_{13}
\left(
  \begin{matrix}
  1 & 0 & -1 & 0 \\
  0 & 0 &  0 & 0 \\
 -1 & 0 &  1 & 0 \\
  0 & 0 &  0 & 0
\right), \qquad s_{13}=\frac{2.0e11 \times 0.005}{2}
```

TO FILL...
