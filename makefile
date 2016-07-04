all:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp main.cpp -D debug="if(1)" -g -Og -Wpedantic -pedantic-errors -Wall -Wextra -Wno-unused-parameter -Werror -std=c++11 2>err.log -o prog || nano err.log
	
clean:
	rm prog

fast:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp main.cpp -Ofast -Wall -Werror -std=c++11 2>err.log -o prog || nano err.log

gprof:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp main.cpp -Ofast -pg -Wall -Werror -std=c++11 -o prog 2>err.log || nano err.log
