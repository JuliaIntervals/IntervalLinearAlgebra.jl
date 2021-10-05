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
\right)
```
