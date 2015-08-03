# Laguerres-Method
This is a suite of Octave functions for graphing the basin of convergence under Laguerre's Method by
calling BasinLaguerreWithCycles.m.  Myriad things will be output, such as the roots and 
nonconvergent cycles.  The resulting graph will be an Argand diagram, with each pixel being an initial 
approximation to the root for Laguerre's Method, and with the pixel being colored according to 
what it eventually converges to.  The hot() colormap, using yellows, oranges, reds, and deep browns, 
correspond to the point eventually converging to a root, with the roots themselves marked 
as black dots.  The ocean() colormap, using blues, corresponds to the point eventually converging to a two-cycle, 
with the points of the two-cycle being marked with + symbols.  The grey() colormap corresponds to the point
eventually converging to a three-cycle, with the points of the three cycle being marked with triangles.  The copper()
colormap corresponds to the point converging to a four cycle, with the points of the four-cycle being marked as squares.
