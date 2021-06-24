# Learn how to solve an interval linear system using Oettli-Präger theorem by drawing the package logo

using IntervalArithmetic, IntervalLinearAlgebra, IntervalConstraintProgramming, Plots

vars = (:x, :y) # also works using ModelingToolkit @variables
A = [2..4 -2..1; -1..2 2..4]
b = [-2..2, -2..2]

X = IntervalBox(-5..5, 2)

α = 0.8
αo = 1
lw = 0.01
lwo = 1

A1 = [2..4 0..1; 0..2 2..4]
b1 = [0..2, -2..0]

p1 = oettli(A1, b1, X, vars)

plot(p1.inner, legend=false, ratio=1, lw=lw, axis=nothing, 
            border=:none, color="#9558B2", alpha=α, background_color=:transparent)

plot!(p1.boundary, color="#9558B2", lw = lwo, lc="#9558B2", α=αo)
A2 = [2..4 -2..0; -1..0 2..4]
b2 = [0..2, 0..2]

p2 = oettli(A2, b2, X, vars)

plot!(p2.inner, legend=false, lw=lw, color="#389826", alpha=α)
plot!(p2.boundary, color="#389826", lw = lwo, lc="#389826", α=αo)

A3 = [2..4 -2..0; -1..0 2..4]
b3 = [-2..0, -2..0]

p3 = oettli(A3, b3, X, vars)

plot!(p3.inner, legend=false, lw=lw, color="#CB3C33", alpha=α)
plot!(p3.boundary, color="#CB3C33", lw = lwo, lc="#CB3C33", α=αo)

A4 = [2..4 0..1; 0..2 2..4]
b4 = [-2..0, 0..2]

p4 = oettli(A4, b4, X, vars)

plot!(p4.inner, legend=false, lw=lw, color="#4063D8", alpha=α)
plot!(p4.boundary, color="#4063D8", lw = lwo, lc="#4063D8", α=αo)

##