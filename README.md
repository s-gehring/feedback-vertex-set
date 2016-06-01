# Compiling

Use the makefile (just type 'make' as a command in shell). Alternatively, type the following command:

g++ fvs_solver.cpp main.cpp -I ./ -O0 -Wall -Werror -std=c++11

The makefile has the advantage, that the project always gets compiled in the same way on every computer. Other than that, the makefile actually compiles and starts the program immediately. All errors (compiling or executing) are written in the file err.log, which gets opened if there were any errors. If there weren't any errors, the normal program output should appear.

