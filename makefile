all:
	clear
	clear
	g++ fvs_solver.cpp main.cpp -O0 -Wall -Werror -std=c++11 2>err.log -o prog && ./prog 2>err.log || nano err.log
	
clean:
	rm prog

fast:
	clear
	clear
	g++ fvs_solver.cpp main.cpp -O3 -Wall -Werror -std=c++11 2>err.log -o prog && ./prog 2>err.log || nano err.log


