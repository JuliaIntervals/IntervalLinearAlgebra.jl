using IntervalArithmetic, StaticArrays, IntervalLinearAlgebra
import IntervalRootFinding:  gauss_seidel_interval, gauss_elimination_interval, gauss_seidel_contractor# not to overload \ in base

A = @SMatrix [4..6 -1..1 -1..1 -1..1;-1..1 -6.. -4 -1..1 -1..1;-1..1 -1..1 9..11 -1..1;-1..1 -1..1 -1..1 -11.. -9]

b = @SVector [-2..4, 1..8, -4..10, 2..12]

jac = Jacobi()
gs = GaussSeidel()
hbr = HansenBliekRohn()

@btime solve($A, $b, $gs)

@btime solve($A, $b, $jac)

@btime solve($A, $b, $hbr)

# comparison with IntervalRootFinding and base

@btime gauss_seidel_interval($A, $b)
@btime gauss_seidel_contractor($A, $b) #NOTE: THIS IS THE JACOBI method
@btime gauss_elimination_interval($A, $b) # this is the one IRF.jl uses to overload \

@btime $A\$b