# # A simple FEM example
# In this section, a problem based on Example 4.1 from [https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/files/7271616/skalna2006.pdf] is considered. The matrix of the system is obtained using the Finite Element Method.
#
# \fig{lgbm_hp1.svg}
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
# where $s$ is the stiffness, given by $EA/L$, with $E$ being the young modulus, $A$ the area of the cross-section and $L$ the length of that truss element.
#
# The change-of-basis matrix is given by
# ```math
# _G(Q)_L = Q =
# \left(
#   \begin{matrix}
# cos(alpha) & -sin(alpha) & 0 & 0 \\
#  sin(alpha) & cos(alpha) &  0 & 0 \\
#  0 & 0 &  cos(alpha) & sin(alpha) \\
#  0 & 0 &  sin(alpha) & cos(alpha)
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
# For element 1->3 the global stiffness matrix is
# ```math
# K_{13} = s_{13}
# \left(
#   \begin{matrix}
#   1 & 0 & -1 & 0 \\
#   0 & 0 &  0 & 0 \\
# -1 & 0 &  1 & 0 \\
#   0 & 0 &  0 & 0
# \right), \qquad s_{13}=\frac{2.0e11 \times 0.005}{2}
# ```
#
# For element 1->2 the global stiffness matrix is
# ```math
# K_{12} = s_{12} \frac{1}{2}
# \left(
#   \begin{matrix}
#   1 & -1 & -1 & 1 \\
#   -1 & 1 &  1 & -1 \\
# -1 & 1 &  1 & -1 \\
#   1 & -1 &  -1 & 1
# \right), \qquad s_{12}=\frac{2.0e11 \times 0.005}{2 \sqrt{2}}
# ```
#
# TO FILL...

# Parameters
E = 2e11 ; # Young modulus
A = 5e-3 ; # Cross-section area

nodesCMatrix = [ 0. 0. ;
                 1. 1. ;
                 2. 0. ;
                 3. 1. ;
                 4. 0. ]

connecMatrix = [ 1 2 ;
                 1 3 ;
                 2 3 ;
                 2 4 ;
                 3 4 ;
                 3 5 ;
                 4 5 ]
fixedDofs     = [2 9 10 ]

# Functions
# computes the unitary stiffness matrix for an element defined by the coordinates of its first and second nodes
function unitaryStiffnessMatrix( coordFirstNode, coordSecondNode  )
  diff      = (coordSecondNode - coordFirstNode)
  length   = sqrt( diff'*diff )
  c        = diff[1] / length ;
  s        = diff[2] / length ;
  Qloc2glo = [ c -s 0 0 ; s c 0 0 ; 0 0 c -s ; 0 0 s c ] ;
  Kloc     = [ 1 0 -1 0 ; 0 0 0 0 ; -1 0 1 0 ; 0 0 0 0 ] ;
  Kglo     = Qloc2glo * Kloc * transpose(Qloc2glo)
  return     Kglo, length
end

numNodes = size( nodesCMatrix )[1]
numElems = size( connecMatrix )[1]
freeDofs = zeros(Int8, 2*numNodes-length(fixedDofs))
indDof  = 1 ;
counter = 0 ;
while indDof <= (2*numNodes)
  if !(indDof in fixedDofs)
    global counter = counter + 1 ;
    freeDofs[ counter ] = indDof ;
    print(indDof)
  end
  global indDof = indDof + 1 ;
end
print(freeDofs)
KG = zeros( 2*numNodes, 2*numNodes ) ;
FG = zeros( 2*numNodes )

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
KG = KG[ freeDofs, : ]
KG = KG[ :, freeDofs ]
FG = FG[ freeDofs ]

u = KG \ FG

print(u)
