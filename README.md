# Compiling

Use the makefile (just type 'make' as a command in shell). Alternatively, type the following command:

g++ fvs_solver.cpp main.cpp -I ./ -O0 -Wall -Werror -std=c++11

The makefile has the advantage, that the project always gets compiled in the same way on every computer. Other than that, the makefile actually compiles and starts the program immediately. All errors (compiling or executing) are written in the file err.log, which gets opened if there were any errors. If there weren't any errors, the normal program output should appear.

# NewGraph

I finished working on the branch NewGraph and merged everything into master. The algorithm should work now. Boost is still necessary, because the graph structure needs flat_set, hash and property_map... or needed that at some point in development and were forgotten. Who cares?
