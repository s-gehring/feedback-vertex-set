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
#include "graph.hpp"

using namespace std;
using namespace FvsGraph;

namespace fvs {

	

	bool has_cycle(Graph& g);
	
	pair<list<Node>, bool> find_semidisjoint_cycle(Graph& g);
	bool edge_exists_between(Graph& g, Node u, Node v);

	Node get_lowest_degree_node(Graph& g, const set<Node>& u);
	bool creates_circle(Graph& g, const set<Node>& u, const Node& v);
	Node two_neighbour_node(Graph& g, const set<Node> &u, const set<Node>& v); //Removed second set. Some1 previously forgot that?

	void cleanup(Graph& g);
	set<Node> two_approx_fvs(Graph& orig);
	bool is_fvs(const Graph& g, const set<Node>& fvs);

	pair<set<Node>, bool> compute_fvs(Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k); ///< DESCRIPTIVE NAMES PLS

	void read_graph(Graph& g, const char* filepath);
	void print_graph(Graph& g);
  
	void induced_subgraph(Graph& s, Graph& g, const set<Node>& u);
	
	enum direction_tag { forward, inverse };
	void maintain_integrity(Graph& g, set<Node>& u, Node aDeletedNode, direction_tag dt = forward);

}

#endif /* FVS_SOLVER_H */
