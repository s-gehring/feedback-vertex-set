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
    ./prog graphs/biggest_example.graph
    gprof prog > gprof_output.txt
    nano gprof_output.txt


# NewGraph

I finished working on the branch NewGraph and merged everything into master. The algorithm should ~~work now~~ be fixed to actually work. Boost is still necessary, because the graph structure needs flat_set, hash and property_map... or needed that at some point in development and were forgotten. Who cares?
