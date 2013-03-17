PDE Solvers for the Laplace's Equation
--------------------------------------

Description:
	A collection of 4 programs to calculate partial differential equations 		also known as PDE Solvers.

Switching between programs:
	The directory includes a git folder for version control. 
	Each program is represented by using a sepparate branch inside git.
	
	The branches are:
		master ( jacobi)
		parallel_pde_solver (parallel jacobi)
		multi-grid 
		parallel-multi-grid

	The following command switches branch, make sure that git is installed.
	
		git checkout <branch name> 

Build:
	make release 
	make debug ( this is default )
		
Run:
	./PDE_SOLVER	<grid size> <number of iterations> <number of workers> <output file>
	./PARALLEL_PDE_SOLVER <grid size> <number of iterations> <number of workers> <output file>
	./MULTI_GRID <coarsest grid size> <number of iterations> <number of workers> <output file>
	./PARALLEL_MULTI_GRID <coarsest grid size> <number of iterations> <number of workers> <output file>

Requirements:
	Windowns: msvc11 ( VS2012 )
	Linux/unix: GCC 4.4.7

