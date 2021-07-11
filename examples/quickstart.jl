using IntervalLinearAlgebra, IntervalConstraintProgramming, Plots

A = [2..4 -1..1;-1..1 2..4]
b = [-2..2, -1..1]

Xenclose = solve(A, b)
Xexact = solve(A, b, NonLinearOettliPrager())

plot(Xexact.inner, ratio=1, label="exact", legend=:top)
plot!(IntervalBox(Xenclose), label="enclosure")
