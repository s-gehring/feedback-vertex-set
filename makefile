all:
	clear
	clear
	g++ graph.cpp debugger.cpp fvs_solver.cpp main.cpp -I ./ -D __DEBUG -g -Og -Wpedantic -pedantic-errors -Wall -Wextra -Werror -std=c++11 2>err.log -o prog && ./prog 2>err.log || nano err.log
	
clean:
	rm prog

fast:
	clear
	clear
	g++ graph.cpp debugger.cpp fvs_solver.cpp main.cpp -I ./ -Ofast -Wall -Werror -std=c++11 2>err.log -o prog && ./prog 2>err.log || nano err.log

gprof:
	clear
	clear
	g++ graph.cpp debugger.cpp fvs_solver.cpp main.cpp -I ./ -Ofast -pg -Wall -Werror -std=c++11 2>err.log || nano err.log


