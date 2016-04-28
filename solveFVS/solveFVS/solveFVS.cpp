// solveFVS.cpp : Defines the entry point for the console application.
//

#include "stdafx.h" // Dunno lol.
#include <iostream>
#include <boost/graph/graph_traits.hpp> // For defining what a node/edge is.
#include <boost/graph/adjacency_list.hpp> // For defining graphs.
#include <fstream> // Only for reading from file, which is not needed in the finished program.

using namespace std;

/*
** A Graph is basically an adjacency list, which saves connections between
** vertices (vecS) and other vertices (vecS) and in this case is undirected.
** It's not too important to choose undirected in this case, but it might be 
** a bit faster.
*/
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;  


Graph* readGraph(char* filename) {
    // This function must be rewritten, since we expect our graph
    // from the stdio, not a file. Furthermore, each vertex can 
    // consist of numbers, characters and even underlines.
    ifstream infile(filename);
    int a, b;
    Vertex u, v;
    Graph G = Graph(2);
    while(infile >> a >> b) {
        u = vertex(a, G);
        v = vertex(b, G);
        add_edge(u, v, G);

    }
    
    return &G;
}

int main(int argc, char** argv)
{
    // We use the first instance if no other is given via command line argument.
    char* filename = "C:\\Users\\Simon\\Downloads\\pace16-fvs-instances-20160301.tar\\pace16-fvs-instances\\001.graph";
    if(argc>1) { filename = argv[1]; }

    Graph* G = readGraph(filename);
    
    

    cin.ignore(cin.rdbuf()->in_avail());
    cin.get();
    return 0;
}
