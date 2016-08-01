all:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -D debug="if(1)" -g -O0 -Wpedantic -pedantic-errors -Wall -Wextra -Wno-unused-parameter -Werror -D_GLIBCXX_DEBUG -std=c++11 2>err.log -o prog || nano err.log

fast:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -std=c++11 2>err.log -o prog -Ofast -Wall -Werror || nano err.log

	
clean:
	rm prog

gprof:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -Ofast -pg -Wall -Werror -std=c++11 -o prog 2>err.log || nano err.log
