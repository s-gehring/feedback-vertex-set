
#ifndef FVS_SOLVER_H
#define FVS_SOLVER_H

#include <iostream>
#include <set>
#include <vector>
#include <stack>
#include "graph.hpp"

using namespace std;
using namespace FvsGraph;

namespace fvs {
	/**
	* @brief Checks, wether the given graph contains a cycle.
	*
	* This function uses dfs and some shortcuts to check if the given graph contains cycles.
	*
	* @param [in] g The graph.
	* @returns True, if there is a cycle in g.
	*/
	bool has_cycle(Graph& g);

	/**
	* @brief Checks, whether a given graph contains a semidisjoint cycle.
	*
	* This function checks if a given graph contains a semidisjoint cycle,
	* i.e. every vertex of the cycle has degree 2 with at most one exeption.
	*
	* We need this for the 2-approximation algo.
	* @param [in] g The graph.
	* @returns The set of nodes forming the semidisjoint cycle in g and an indicator for the existence.
	*/
	pair<list<Node>, bool> find_semidisjoint_cycle(Graph& g);

	/**
	* @brief Finds the lowest degree node in the graph induced subgraph of a set of nodes.
	*
	* This function iterates over all nodes to find the one with lowest degreen in U.
	* Only neighbours also in u are counted.
	*
	* During the algorithm, it is assumed, that this method always returns a node with
	* at most degree 1.
	*
	* @param [in] g The graph.
	* @param [in] u The set of nodes the induced subgraph is created from.
	* @returns The lowest degree node in u along with its degree.
	*/
	Node get_lowest_degree_node(const Graph& g, const set<Node>& u);


	/**
	* @brief This function checks, if two of the neighbours of v belong to the same tree in g[u].
	*
	* This check is done by checking wether any neighbour of v in u has a neighbour (also in u), that
	* is also neighbour to v. If this happens, it is assumed that they belong in the same tree, thus
	* the function returns true.
	*
	* @param [in] g The basic graph.
	* @param [in] u The nodeset that does not induce the subgraph for which were checking the neighbourhood of v.
	* @param [in] v A node, which might connect a circle in g[u].
	* @returns True, if a neighbour of a neighbour of v is a neighbour of v.
	*/
	bool creates_circle(Graph& g, const set<Node>& u, const Node& v);

	/**
	* @brief Finds a node in u, which has atleast two neighbours in v with respect to g.
	*
	* This function iterates over all elements of u to find a node which has two neighbours
	* in v, where both u and v are sets of nodes of the same graph g. Returns a node, if one
	* was found or a null_vertex() if no fitting node was found.
	*
	* @param [in] g The graph this is based on.
	* @param [in] u A node subset of g which partitions g in u and g-u.
	* @returns A node of u with atleast two neighbours in v, if such a node exists or a null_vertex()
	*		otherwise.
	*/
	Node two_neighbour_node(Graph& g, const set<Node> &u, const set<Node>& v);

	/**
	* @brief: Deletes all vertices of degree at most 1 along with all incident edges from a given graph.
	*
	* As long as a given graph has vertices of degree at most 1, all incident edges and all vertices
	* are deleted. We need this as a subroutine for the 2-approx-algo.
	*
	* @param [in] g The graph.
	*/
	void cleanup(Graph& g);

	/**
	*@brief: Computes a 2 - approximation of an fvs for a given graph.
	*
	* This function computes a 2 - approximation of a feedback vertex set following
	* Bafna et al, A 2 - APPROXIMATION ALGORITHM FOR THE UNDIRECTED
	* FEEDBACK VERTEX SET PROBLEM, 1999. It runs in O(n ^ 2) if implemented correctly :D.
	*
	* We need this to initialize the fvs for the iterative compression process.
	*
	* @param[in] orig The graph.
	* @returns The feedback vertex set.
	*/
	set<Node> two_approx_fvs(Graph& orig);

	/**
	*@brief: Checks whether a given set is a feedback vertex set in a given graph.
	*
	* @param[in] g The graph.
	* @param[in] fvs The set to be checked.
	* @returns True, if the given set is a feedback vertex set in g.
	*/
	bool is_fvs(const Graph& g, const set<Node>& fvs);

	/**
	*@brief: Computes the union of two given sets.
	*
	* @param[in] S The first set.
	* @param[in] T The second set.
	* @returns The union of the sets.
	*/
	set<Node> set_union(const set<Node> S, const set<Node> T);

	/**
	* @brief Computes a feedback vertex set.
	*
	* This algorithm was proposed in the paper, that we have here:
	*	https://www.dropbox.com/sh/ar26siyo2cjw6y1/AACqJWkA0YHXxkg5FTz2ZEeBa/ChenFLLV08_ImprovedAlgorithmsForFeedbackVertexSetProblems.pdf?dl=0
	*
	* This function right now resembles a direct copy of the pseudo code of the article.
	*
	* @param [in] orig The original graph, that execution started with.
	* @param [in] g The graph to find a feedback vertex set for.
	* @param [in] f A subgraph of g, for which f is a forest and g-f is a forest.
	* @param [in] v2 The second set of the forest bipartition (f, v2).
	* @param [in] k The currently guessed size of a min. feedback vertex set.
	* @returns A pair of a set of nodes and a bool. The set of nodes contains a part of the feedback
	*		vertex set. The bool will be false, if the algorithm decides that there is no fvs.
	*/
	pair<set<Node>, bool> forest_bipartition_fvs(Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k);

	/**
	*@brief: Tries to decrease the size of a given fvs in a given graph by 1.
	*
	* @param[in] orig The graph.
	* @param[in] S The feedback vertex set which will be compressed.
	* @returns A new fvs with size of the old fvs -1 if it exists, otherwise false.
	*/
	pair<set<Node>, bool> compression_fvs(const Graph& orig, const set<Node>& S);

	/**
	*@brief: Computes the minimum feedback vertex set for a given graph.
	*
	* Following Chen et al, Improved Algorithms for FVS problems, this algorithm uses
	* iterative compression to compute a minimum feedback vertex set for a given graph.
	*
	* @param[in] orig The graph.
	* @returns A minimum feedback vertex set.
	*/
	set<Node> compute_min_fvs(const Graph& orig);

	/**
	*@brief: Computes the minimum feedback vertex set for a given graph using brute force.
	*
	* @param[in] orig The graph.
	* @returns A minimum feedback vertex set.
	*/
	set<Node> brute_force_fvs(const Graph& orig);

	/**
	*@brief: Computes the difference set of two given sets.
	*
	* @param[in] S The first set.
	* @param[in] T The set the will be substracted from the first one.
	* @returns The difference of the sets.
	*/
	set<Node> set_minus(const set<Node> S, const set<Node> T);

	/**
	* @brief Reads in a graph from a standard text format.
	*
	* This function reads a graph from a .txt file with the following structure. The first
	* line contains the number of nodes, the second on the number of edges in the graph. All
	* other lines contain to positive integers denoting source and target of an edge.
	*
	* I know that this isnt the format of the challenge, but this will suffice for testing.
	*
	* @param [out] g The graph will be written to this object. All data will be cleared before reading the new graph.
	* @param [in] filepath The path to the file containing the graph.
	*/
	void read_graph(Graph& g, const char* filepath);

	/**
	* @brief Print g to stdout.
	*/
	void print_graph(Graph& g);

	/**
	* Prints a set of nodes to stdout.
	*
	* @param [in] s Set of nodes to be printed.
	*/
	void print_nodes(set<Node>& s);

	/**
	* @brief Creates the induced subgraph g[u].
	*/
	void induced_subgraph(Graph& s, Graph& g, const set<Node>& u);
	
	enum direction_tag { forward, inverse };
}

#endif /* FVS_SOLVER_H */
