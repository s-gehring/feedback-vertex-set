all:
	clear
	clear
	g++  fvs_solver.cpp main.cpp -I ./ -g -Og -Wall -Werror -std=c++11 2>err.log -o prog && ./prog 2>err.log || nano err.log
	
clean:
	rm prog

fast:
	clear
	clear
	g++ fvs_solver.cpp main.cpp -I ./ -Ofast -Wall -Werror -std=c++11 2>err.log -o prog && ./prog 2>err.log || nano err.log

gprof:
	clear
	clear
	g++ fvs_solver.cpp main.cpp -I ./ -Ofast -pg -Wall -Werror -std=c++11 2>err.log || nano err.log


