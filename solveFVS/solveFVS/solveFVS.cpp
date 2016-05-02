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
// TODO:
// Implement all this shit. Preferably Object Oriented.

// Returns true iff G has at least one cycle. In O(|G|) or better.
bool hasCycle(Graph G);                                     

// ... what do you think?
bool hasEdge(Graph G, Vertex u, Vertex v);                  

// For any pair of neighbours (x,y) of u in U: Is there a path in G from x to y,
// only visiting vertices of U?
// This method is going to be called a lot. Better make it fast.
bool sameTreeNeighbours(Graph G, list<Vertex> U, Vertex u); 

// Obviously not making a copy.
Graph removeVertex(Graph G, Vertex x);

// Get a node with a degree of <= 1 in U.
// That means the returned node has at most 1 neighbor in U, but maybe more in G.
// Such a vertex is guaranteed to exist.
Vertex getLowDegreeNode(Graph G, list<Vertex> U);

int getDegree(Graph G, Vertex x);
// We need this!
Graph operator -(Graph G, Vertex x) { // Does operator overloading work this way?
    return removeVertex(G, x);
}
list<Vertex> operator -(list<Vertex> V, Vertex x) {
    V.remove(x);
    return V;
}
list<Vertex> operator +(list<Vertex> V, Vertex y) {
    V.insert(y);
    return V;
}


list<Vertex> Feedback(Graph G, list<Vertex> V_1, list<Vertex> V_2, int k) {
    /*
        https://www.dropbox.com/sh/ar26siyo2cjw6y1/AACqJWkA0YHXxkg5FTz2ZEeBa/ChenFLLV08_ImprovedAlgorithmsForFeedbackVertexSetProblems.pdf?dl=0
        ... I'm so sorry for what I've done here.
    */


    list<Vertex> F;
    if(k < 0 || (k==0 && hasCycle(G))) return F; // Return null or false, not F, not empty list. Dunno how to implement this though.
    if(k >= 0 && !hasCycle(G)) return F;
    
    // Find a Node in V_1 with at least two neighbors in V_2
    Vertex w = (Vertex) nullptr; // Fuck the police.
    for(Vertex v : V_1) {
        int neighbours = 0;
        for(Vertex w : V_2) {
            if(hasEdge(G, v, w)) {
                if(++neighbours == 2) {
                    // Found a Vertex w in V_1 which has at least two neighbors in V_2.
                    w = v;
                    break;
                }
            }
        }
        if(w != (Vertex) nullptr) break;
    }

    /*
        Sorry for that shit name 'w'. It's called w in the paper.
    */

    if(w != (Vertex) nullptr) {
/*3.1.*/if(sameTreeNeighbours(G, V_2, w)) {
            
            F = Feedback(G-w, V_1-w, V_2, k-1);
            if(F == nullptr) return nullptr; // Again, check if F==NULL and return NULL. 
            F.insert(w); // help :(
            return F;
/*3.2.*/} else {

            F = Feedback(G-w, V_1-w, V_2, k-1);
            if(F != nullptr) return F; // If F!=NULL return F;
            return Feedback(G, V_1-w, V_2+w, k);
        }
    } else {
        w = getLowDegreeNode(G, V_1);
/*4.1.*/if(getDegree(G, w) <= 1) {
    return Feedback(G-w, V_1-w, V_2, k); // VS indents this line weirdly.
/*4.2.*/} else {
    return Feedback(G, V_1-w, V_2+w, k);
        }
    }
}



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
