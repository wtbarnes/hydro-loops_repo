hydro-loops
===========
## About
The hydro-loops code computes hydrostatic solutions (i.e. _dt_=0) for coronal loops, 
semi-stable (on the order of several hours) structures in the upper solar atmosphere 
consisting of heated plasma along arcing magnetic field lines rooted in the solar 
chromosphere. In particular, the code solves the hydrostatic equations with steady flows. 
Details of the calculation can be found in <a 
href="http://adsabs.harvard.edu/abs/2001ApJ...550.1036A">Aschwanden et al. (2001)</a>.
## Model Details
The code uses a shooting method to satisfy a boundary condition on the conductive heat 
flux _F(s)_, in particular _F(s=L)=0_ for a loop of half-length _L_. As of now, the 
equations are solved using a simple, non-adaptive Euler solver. 
## Compiling and Executing
