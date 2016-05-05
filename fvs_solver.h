/*
* File:   fvs_solver.h
* Author: fabian
*
* Created on 04. Mai 2016, 15:40
*
* Version: 0.1
* This file describes the interface used to solve the fvs problem.
*/

#ifndef FVS_SOLVER_H
#define FVS_SOLVER_H

#include <iostream>
#include <set>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/graph/graph_traits.hpp"

using namespace boost;
using namespace std;

namespace fvs {

	typedef adjacency_list<vecS, vecS, undirectedS> Graph;
	typedef graph_traits<Graph>::vertex_descriptor Node;
	typedef graph_traits<Graph>::edge_descriptor Edge;


	bool has_cycle(const Graph& g);
	bool has_semidisjoint_cycle(const Graph& g);
	bool edge_exists_between(const Graph& g, Node u, Node v);

	Node get_lowest_degree_node(const Graph& g, const set<Node>& u);
	bool creates_circle(const Graph& g, const set<Node>& u, const Node& v);
	Node two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node> &v);

	void cleanup(Graph& g);
	set<Node> two_approx_fvs(Graph& g);

	pair<set<Node>, bool> compute_fvs(Graph& g, set<Node> f, int k); ///< DESCRIPTIVE NAMES PLS

	void read_graph(Graph& g, const char* filepath);

	void induced_subgraph(Graph& s, const Graph& g, const set<Node>& u);

	/**
	* A dfs visitor that helps to finds circles in a graph. Finds a circle, if a back, forward or cross edge is found.
	*/
	struct CycleVisitor : public default_dfs_visitor {
		CycleVisitor(bool &b) {_circle=&b;};
		CycleVisitor(const CycleVisitor& other) : _circle(other._circle) {};
		void back_edge(Edge e, Graph g) { *_circle = true; };
		void forward_or_cross_edge(Edge e, Graph g) { *_circle = true; };

		bool *_circle; //True if a circle was found.
	};
}

#endif /* FVS_SOLVER_H */
