# # A simple FEM application
#
# ## Introduction
#
# ### General description
# The Finite Element Method is widely used to solve PDEs in Engineering applications and particularly in Structural Analysis problems [[BAT14]](@ref). The procedure consists discretizing the domain into _elements_ and assembly a system of balance equations. For linear problems, this system can be usually written as
#
# ```math
# K \cdot d = f
# \qquad
# K = \sum_e K_e
# ```
# where $f$ is the vector of external loads, $K_e$ is the element stiffness matrix and $d$ is the vector of unknown displacements.
#
# ### FEM for truss structures
# A frequent and simple type of structures are _Truss structures_, which are formed by bars connected but not welded. Truss models are considered during the conceptual design of bridges or other structures.
#
# The stiffness matrix of a truss element in the local coordinate system is given by
# ```math
# K_L = s
# \left(
#  \begin{matrix}
#  1 & 0 & -1 & 0 \\
#  0 & 0 &  0 & 0 \\
# -1 & 0 &  1 & 0 \\
#  0 & 0 &  0 & 0
# \end{matrix}
# \right),
# ```
#
# where $s =\frac{E A}{L}$, with $E$ being the Young modulus, $A$ the area of the cross-section and $L$ the length of that truss element.
#
# The change-of-basis matrix is given by
# ```math
# _G(Q)_L = Q =
# \left(
#   \begin{matrix}
# \cos(\alpha) & -\sin(\alpha) & 0 & 0 \\
#  \sin(\alpha) & \cos(\alpha) &  0 & 0 \\
#  0 & 0 &  \cos(\alpha) & -\sin(\alpha) \\
#  0 & 0 &  \sin(\alpha) & \cos(\alpha)
# \end{matrix}
# \right),
# ```
#
# The system of equations for each element is written in local coordinates as
# ```math
# K_L d_L = f_L
# ```
# and using the change-of-basis we obtain
# ```math
# K_G d_G = f_G \qquad K_G = Q K_L Q^T
# ```
#
# The unitary stiffness matrix (for $s=1$) can be computed using the following function.
function unitaryStiffnessMatrix( coordFirstNode, coordSecondNode  )
  diff      = (coordSecondNode - coordFirstNode)
  length   = sqrt( diff' * diff )
  c        = diff[1] / length ;
  s        = diff[2] / length ;
  Qloc2glo = [ c -s 0 0 ; s c 0 0 ; 0 0 c -s ; 0 0 s c ] ;
  Kloc     = [ 1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0 ] ;
  Kglo     = Qloc2glo * Kloc * transpose(Qloc2glo)
  return     Kglo, length
end
#
# ## Example problem
#
# A problem based on Example 4.1 from [[SKA06]](@ref) is considered. The following diagram shows the truss structure considered.
#
# ```@raw html
# <img src="../../assets/trussDiagram.svg" style="width: 100%" alt="truss diagram"/>
# ```
#
# ### Case with fixed parameters
#
# The scalar parameters considered are given by
E = 2e11 ; # Young modulus
A = 5e-3 ; # Cross-section area
# The coordinate matrix is given by
## coordinates   x  y
nodesCMatrix = [ 0. 0. ;
                 1. 1. ;
                 2. 0. ;
                 3. 1. ;
                 4. 0. ];
# the connectivity matrix is given by
## connectivity  start end
connecMatrix = [ 1     2 ;
                 1     3 ;
                 2     3 ;
                 2     4 ;
                 3     4 ;
                 3     5 ;
                 4     5 ];
# and the fixed degrees of freedom (supports) are defined by the vector
fixedDofs     = [2 9 10 ];
#
# calculations
numNodes = size( nodesCMatrix )[1]; # compute the number of nodes
numElems = size( connecMatrix )[1]; # compute the number of elements
freeDofs = zeros(Int8, 2*numNodes-length(fixedDofs));
indDof  = 1 ; counter = 0 ;
while indDof <= (2*numNodes)
  if !(indDof in fixedDofs)
    global counter = counter + 1 ;
    freeDofs[ counter ] = indDof ;
  end
  global indDof = indDof + 1 ;
end
#
# assembly
KG = zeros( 2*numNodes, 2*numNodes );
FG = zeros( 2*numNodes );
for elem in 1:numElems
  print(" assembling stiffness matrix of element ", elem , "\n")
  indexFirstNode  = connecMatrix[ elem, 1 ]
  indexSecondNode = connecMatrix[ elem, 2 ]
  dofsElem = [2*indexFirstNode-1 2*indexFirstNode 2*indexSecondNode-1 2*indexSecondNode ]
  KGelem, lengthElem = unitaryStiffnessMatrix( nodesCMatrix[ indexSecondNode, : ], nodesCMatrix[ indexFirstNode, : ] )
  stiffnessParam = E * A / lengthElem ;
  for i in 1:4
    for j in 1:4
      KG[ dofsElem[i], dofsElem[j] ] = KG[ dofsElem[i], dofsElem[j] ] + stiffnessParam * KGelem[i,j]
    end
  end
end
FG[4] = -1e4 ;
KG = KG[ freeDofs, : ] ;
KG = KG[ :, freeDofs ] ;
FG = FG[ freeDofs ];

u = KG \ FG
UG = zeros( 2*numNodes );
UG[ freeDofs ] = u ;
#
# #### Deformed structure
#
# The reference (dashed blue line) and deformed (solid red)  configurations of the structure are ploted. Since the displacements are very small, a `scaleFactor` is considered to amplify the deformation and ease the visualization.
#
using Plots
scaleFactor = 2e3 ;
plot();
for elem in 1:numElems
  indexFirstNode  = connecMatrix[ elem, 1 ];
  indexSecondNode = connecMatrix[ elem, 2 ];
  ## plot reference element
  plot!( nodesCMatrix[ [indexFirstNode, indexSecondNode], 1 ],
         nodesCMatrix[ [indexFirstNode, indexSecondNode], 2 ],
         linestyle = :dash,  aspect_ratio = :equal,
         linecolor = "blue", legend = false)

  ## plot deformed element
  plot!( nodesCMatrix[ [indexFirstNode, indexSecondNode], 1 ]
           + scaleFactor* [ UG[indexFirstNode*2-1], UG[indexSecondNode*2-1]] ,
         nodesCMatrix[ [indexFirstNode, indexSecondNode], 2 ]
           + scaleFactor* [ UG[indexFirstNode*2  ], UG[indexSecondNode*2  ]] , markershape = :circle, aspect_ratio = :equal, linecolor = "red",
           linewidth=1.5, legend = false )
end
xlabel!("x (m)") # hide
ylabel!("y (m)") # hide
title!( "Deformed with scale factor " * string(scaleFactor) ) # hide
savefig("deformed.png") # hide
#
# ![](deformed.png)
#
# ### Problem with interval parameters

md"""
Suppose now we have a 10% uncertainty for the stiffness $s_{23}$ associated with the third
element. To model the problem, we introduce the symbolic variable `s23` using the IntervalLinearAlgebra macro `@linvars`.
"""

using IntervalLinearAlgebra
@linvars s23

# now we can construct the matrix as before

KGp = zeros(AffineExpression{Float64}, 2*numNodes, 2*numNodes );
for elem in 1:numElems
  print(" assembling stiffness matrix of element ", elem , "\n")
  indexFirstNode  = connecMatrix[ elem, 1 ]
  indexSecondNode = connecMatrix[ elem, 2 ]
  dofsElem = [2*indexFirstNode-1 2*indexFirstNode 2*indexSecondNode-1 2*indexSecondNode ]
  KGelem, lengthElem = unitaryStiffnessMatrix( nodesCMatrix[ indexSecondNode, : ], nodesCMatrix[ indexFirstNode, : ] )
  if elem == 3
    stiffnessParam = s23
  else
    stiffnessParam = E * A / lengthElem ;
  end
  for i in 1:4
    for j in 1:4
      KGp[ dofsElem[i], dofsElem[j] ] = KGp[ dofsElem[i], dofsElem[j] ] + stiffnessParam * KGelem[i,j]
    end
  end
end
KGp = KGp[ freeDofs, : ]
KGp = KGp[ :, freeDofs ]

# Now we can construct the [`AffineParametricArray`](@ref)

KGp = AffineParametricArray(KGp)

# The range of the stiffness is

srange = E * A / sqrt(2) ± 0.1 * E * A / sqrt(2)

# To solve the system, we could of course just subsitute `srange` into the parametric matrix
# `KGp` and solve the "normal" interval linear system

usimple = solve(KGp(srange), Interval.(FG))

# This approach, however soffers from the [dependency problem](https://en.wikipedia.org/wiki/Interval_arithmetic#Dependency_problem)
# and hence the compute displacement will be an overestimation of the true displacement.

# To mitigate this issue, algorithms to solve linear systems with parameters have been developed.
# In this case we use the algorithm presented in [[SKA06]](@ref)

uparam = solve(KGp, FG, srange)

# We can now compare the naive and parametric solution

hcat(usimple, uparam)/1e-6
