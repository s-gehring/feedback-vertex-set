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

	/**
	* A dfs visitor that helps to find semi-disjoint cycles in a graph.
	*/
	struct SemiDisjointCycleVisitor : public dfs_visitor<> {
		SemiDisjointCycleVisitor(bool &b, set<Node> &sdc) : _sdcycle(b), _sdVertices(sdc) {};
		SemiDisjointCycleVisitor(const SemiDisjointCycleVisitor& other) : _sdcycle(other._sdcycle), _sdVertices(other._sdVertices) {};
		void tree_edge(Edge e, const Graph& g) {
			//			cout << "Found edge " << source(e,g) << " -> " << target(e,g) << endl;
			foundEdges.push_back(make_pair(source(e, g), target(e, g)));
		};
		void back_edge(Edge e, const Graph& g) {
			//			cout << "Back edge: " << source(e, g) << " -> " << target(e,g) << endl;

			if (find(foundEdges.begin(), foundEdges.end(), make_pair(target(e, g), source(e, g))) == foundEdges.end()) {
				//				cout << "Not yet found." << endl;
				foundEdges.push_back(make_pair(source(e, g), target(e, g)));
				// found cycle, check if it is semi-disjoint
				if (!_sdcycle) {
					// get vertices of cycle going backwards along the edges
					vector<pair<Node, Node>>::iterator it = foundEdges.end()-1;
					int vertices_with_high_degree = 0;
					while((*it).first != target(e, g)) {
						_sdVertices.insert((*it).second);
						// count number of vertices with degrees >2
						if (in_degree((*it).second, g) > 2) {
							vertices_with_high_degree++;
						}
						--it;
					}
					_sdVertices.insert((*it).second);
					if (in_degree((*it).second, g) > 2) {
						vertices_with_high_degree++;
					}
					if(vertices_with_high_degree <= 1) {
						_sdcycle = true;
					}
					else{
						_sdVertices.clear();
					}
				}
			}
		};

		bool& _sdcycle; //True if a semi-disjoint-cycle was found.
		vector<pair<Node, Node>> foundEdges;
		set<Node>& _sdVertices;
	};
}

#endif /* FVS_SOLVER_H */
