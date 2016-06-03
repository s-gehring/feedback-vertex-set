/*
* File:   fvs_solver.h
*
* This file describes the interface used to solve the fvs problem.
*/

#ifndef FVS_SOLVER_H
#define FVS_SOLVER_H

#include <iostream>
#include <set>
#include <vector>
#include <stack>
/*
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/undirected_dfs.hpp"
#include "boost/graph/graph_traits.hpp"
*/
#include "graph.cpp"

using namespace boost;
using namespace std;
using namespace Graph;

namespace fvs {

	

	bool has_cycle(const Graph& g);
	pair<set<Node>, bool> find_semidisjoint_cycle(const Graph& g);
	bool edge_exists_between(const Graph& g, Node u, Node v);

	Node get_lowest_degree_node(const Graph& g, const set<Node>& u);
	bool creates_circle(const Graph& g, const set<Node>& u, const Node& v);
	Node two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node>& v); //Removed second set. Some1 previously forgot that?

	void cleanup(Graph& g, map<Node, double>& weights);
	set<Node> two_approx_fvs(const Graph& orig);

	pair<set<Node>, bool> compute_fvs(const Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k); ///< DESCRIPTIVE NAMES PLS

	void read_graph(Graph& g, const char* filepath);
	void print_graph(const Graph& g);

	void induced_subgraph(Graph& s, const Graph& g, const set<Node>& u);
	
	enum direction_tag { forward, inverse };
	void maintain_integrity(Graph& g, set<Node>& u, Node aDeletedNode, direction_tag dt = forward);

	
}

#endif /* FVS_SOLVER_H */
