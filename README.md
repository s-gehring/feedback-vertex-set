# Compiling

Use the makefile. 

    make          
Compile for debugging. No speed ups, include debug symbols for gdb.

    make fast     
Compile for speed tests only. This is basically how our program will be compiled at the end.

    make gprof    
Compile for speed profiling with gprof (or similar). Uses flag -pg and debug symbols.

# gprof

    sudo apt-get install gprof
    make gprof
    ./prog < graphs/biggest_example.graph
    gprof prog > gprof_output.txt
    nano gprof_output.txt


# Program usage

    ./prog < [path_to_graph_file]
