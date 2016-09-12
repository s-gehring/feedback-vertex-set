default:
	g++ graph.cpp fvs_solver.cpp degree3.cpp galois.cpp gauss.cpp lin_parity.cpp matrix.cpp matrix2.cpp static.cpp bin_count.cpp main.cpp -std=c++11 2>err.log -o prog -Ofast -Wall -Werror -pedantic || nano err.log
clean:
	rm prog
