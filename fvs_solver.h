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
#include "boost/graph/undirected_dfs.hpp"
#include "boost/graph/graph_traits.hpp"

using namespace boost;
using namespace std;

namespace fvs {

	typedef adjacency_list<vecS, vecS, undirectedS,
		no_property,
		property<edge_color_t, default_color_type>> Graph;
	typedef graph_traits<Graph>::vertex_descriptor Node;
	typedef graph_traits<Graph>::edge_descriptor Edge;


	bool has_cycle(const Graph& g);
	bool has_semidisjoint_cycle(const Graph& g);
	bool edge_exists_between(const Graph& g, Node u, Node v);

	Node get_lowest_degree_node(const Graph& g, const set<Node>& u);
	bool creates_circle(const Graph& g, const set<Node>& u, const Node& v);
	Node two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node>& v); //Removed second set. Some1 previously forgot that?

	void cleanup(Graph& g);
	set<Node> two_approx_fvs(Graph& g);

	pair<set<Node>, bool> compute_fvs(const Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k); ///< DESCRIPTIVE NAMES PLS

	void read_graph(Graph& g, const char* filepath);
	void print_graph(const Graph& g);

	void induced_subgraph(Graph& s, const Graph& g, const set<Node>& u);
	
	enum direction_tag { forward, inverse };
	void maintain_integrity(Graph& g, set<Node>& u, Node aDeletedNode, direction_tag dt = forward);

	/**
	* A dfs visitor that helps to finds circles in a graph. Finds a circle, if a back, forward or cross edge is found.
	*/
	struct CycleVisitor : public dfs_visitor<> {
		CycleVisitor(bool &b) : _circle(b) { };
		CycleVisitor(const CycleVisitor& other) : _circle(other._circle) {};
		void tree_edge(Edge e, const Graph& g) {
//			cout << "Found edge " << source(e,g) << " -> " << target(e,g) << endl;
			foundEdges.insert(make_pair(source(e,g), target(e,g)));
		}
		void back_edge(Edge e, const Graph& g) {
//			cout << "Back edge: " << source(e, g) << " -> " << target(e,g) << endl;

			if (foundEdges.end() == foundEdges.find(make_pair(target(e,g), source(e,g)))) {
//				cout << "Not yet found." << endl;
				_circle = true;
				foundEdges.insert(make_pair(source(e,g), target(e,g)));
			}
		};

		bool& _circle; //True if a circle was found.
		set<pair<Node, Node>> foundEdges;
	};
}

#endif /* FVS_SOLVER_H */
