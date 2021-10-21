# # A simple FEM application
#
# ## Introduction
# The Finite Element Method is widely used to solve PDEs in Engineering applications and particularly in Structural Analysis problems [Add citation Bathe]. A specific case os structures is Truss structures, where trusses or bars are connected but not welded. Truss models are considered during the conceptual design of bridges or other structures.
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
# where $s$ is the stiffness, given by $EA/L$, with $E$ being the Young modulus, $A$ the area of the cross-section and $L$ the length of that truss element.
#
# The change-of-basis matrix is given by
# ```math
# _G(Q)_L = Q =
# \left(
#   \begin{matrix}
# \cos(\alpha) & -\sin(\alpha) & 0 & 0 \\
#  sin(alpha) & cos(alpha) &  0 & 0 \\
#  0 & 0 &  cos(alpha) & sin(alpha) \\
#  0 & 0 &  sin(alpha) & cos(alpha)
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
# The stiffness matrix can be computed using the following function, which computes the unitary stiffness matrix for an element defined by the coordinates of its first and second nodes.
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
#
# ## Problem with fixed parameters
# In this section, a problem based on Example 4.1 from [https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl/files/7271616/skalna2006.pdf] is considered. The following diagram shows the truss structure considered.
#
# \fig{../../assets/trussDiagram.png}
#
# The scalar parameters considered are given by
E = 2e11 ; # Young modulus
A = 5e-3 ; # Cross-section area
# while the coordinate matrix is given by
nodesCMatrix = [ 0. 0. ;
                 1. 1. ;
                 2. 0. ;
                 3. 1. ;
                 4. 0. ]
# and connectivity matrix is given by
connecMatrix = [ 1 2 ;
                 1 3 ;
                 2 3 ;
                 2 4 ;
                 3 4 ;
                 3 5 ;
                 4 5 ]
# and the fixed degrees of freedom (supports) are
fixedDofs     = [2 9 10 ]

# calculations
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

# assembly
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

# ## Problem with interval parameters
