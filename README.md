#hydro-loops
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
A `makefile` is included to compile the code into the executable `hydro-loops`. Download 
the compressed file or fork the repository (`git clone 
https://github.com/wtbarnes/hydro-loops_repo.git`) and then run `make` in the 
`hydro-loops_repo' directory. 
### Configuring Input Parameters
The input file is a simple text file `hydro-loops_parameters.txt`. The structure is as 
follows:

1. Number of grid cells
2. Heating switch (0=uniform, 1=non-uniform)
3. Radiative-loss switch (0=simple, 1=full)
4. Lower bound for heating
5. Upper bound for heating
6. T0--temperature at loop base
7. n0--density at loop base
8. h0--height of loop base above the solar surface
9. v0--velocity at loop base(v0=0 leads to solving the standard hydrostatic equations)
10. Threshold value for _F(s=L)_ boundary condition (i.e. how close to zero will the model try to get)
11. species--electron or ion fluid

Note that the loop half-length and heating scale-height (important only for non-uniform 
heating) are set via command-line arguments. For example, to run the model for _L_=50 Mm 
and _Sh_=100 Mm, one would run `./hydro-loops 50 100`.
