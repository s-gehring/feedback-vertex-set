all:
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -std=c++11 -o prog -O3

deg3:
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp Tests.cpp bin_count.cpp deg3main.cpp -std=c++11 -o prog -O3

debug:
	clear
	clear
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -D debug="if(1)" -g -Og -Wpedantic -pedantic-errors -Wall -Wextra -Wno-unused-parameter -Werror -std=c++11 2>err.log -o prog || nano err.log

debugi:
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -D debug="if(1)" -D_GLIBCXX_DEBUG -g -O0 -std=c++11 -o prog
