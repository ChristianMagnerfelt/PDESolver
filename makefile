debug: pde_solver.cpp
	g++ -g -O0 -Wall -Wextra -o MULTI_GRID pde_solver.cpp -std=c++0x -fopenmp
	
release: pde_solver.cpp
	g++ -O3 -Wall -Wextra -o MULTI_GRID pde_solver.cpp -std=c++0x -fopenmp
	
clean: 
	rm MULTI_GRID
