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

#include "boost\graph\adjacency_list.hpp"
#include "boost\graph\depth_first_search.hpp"
#include "boost\graph\graph_traits.hpp"

using namespace boost;
using namespace std;

namespace fvs {

	typedef adjacency_list<> Graph;
	typedef graph_traits<Graph>::vertex_descriptor Node;
	typedef graph_traits<Graph>::edge_descriptor Edge;


	bool has_cycle(const Graph& g);
	bool edge_exists_between(const Graph& g, Node u, Node v);

	Node get_lowest_degree_node(const Graph& g, const set<Node>& u);

	set<Node> compute_fvs(Graph& g, set<Node> v1, set<Node> v2, int k); ///< DESCRIPTIVE NAMES PLS

	void read_graph(Graph& g, const char* filepath);

	/**
	* A dfs visitor that helps to finds circles in a graph. Finds a circle, if a back, forward or cross edge is found.
	*/
	struct CycleVisitor : public default_dfs_visitor {
		CycleVisitor() : _circle(false) {};
		CycleVisitor(const CycleVisitor& other) : _circle(other._circle) {};
		void back_edge(Edge e, Graph g) { *_circle = true; };
		void forward_or_cross_edge(Edge e, Graph g) { *_circle = true; };

		bool *_circle; //True if a circle was found.
	};
}

#endif /* FVS_SOLVER_H */