
#include <utility>
#include <limits.h>
#include <fstream>

#include "fvs_solver.h"

using namespace fvs;

/**
* @brief Checks, wether the given graph contains a circle.
*
* This function uses dfs to check if the given graph contains circles. Right now
* a circle is found, when the dfs finds a back, forward or cross edge (since the graph is undirected).
* 
* Boost's visitor concept is used to transport info.
* 
* @param [in] g The graph.
* @returns True, if there is a circle in g (hopefully).
*/
bool fvs::has_cycle(const fvs::Graph& g) {
	CycleVisitor cv;
	depth_first_search(g, visitor(cv));
	return *cv._circle;
}

/**
* This is basically an example on how to do that with boost.
* Dont use this method, just call the original.
*/
bool fvs::edge_exists_between(const Graph& g, Node u, Node v) {
	return edge(u, v, g).second;
}

/**
* @brief Finds the lowest degree node in the graph.
*
* This function iterates over all nodes to find the one with lowest degreen in U.
* Only neighbours also in u are counted.
*
* During the algorithm, it is assumed, that this method always returns a node with 
* degree 1.
*
* PLEASE SOMEBODY COME UP WITH A NAME FOR U. A DESCRIPTIVE ONE!
*
* @param [in] g The graph.
* @param [in] u A set of nodes. 
* @returns The lowest degree node in u.
*/
Node fvs::get_lowest_degree_node(const Graph &g, const set<Node>& u) {
	if (u.size() == 0) {
		throw runtime_error("Error: Searching lowest degree node in an empty set.");
	}

	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
	pair<size_t, Node> lowestDegreeNode = make_pair((size_t) - 1, NULL);

	for (const auto& i : u) {
		if (u.end() != u.find(i)) {
			size_t edgesToOtherUs = 0;
			pair<edge_iterator, edge_iterator> eIt = out_edges(i, g);
			for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
				if (u.end() != u.find(target((*it), g))) {
					++edgesToOtherUs;
				}
			}

			if (edgesToOtherUs < lowestDegreeNode.first) {
				lowestDegreeNode = make_pair(edgesToOtherUs, i);
			}
		}
	}
	return lowestDegreeNode.second;
}

/**
* @brief This function checks, if two of the neighbours of v belong to the same tree in g[u].
*
* This check is done by checking wether any neighbour of v in u has a neighbour (also in u), that
* is also neighbour to v. If this happens, it is assumed that they belong in the same tree, thus
* the function returns true.
*
* @param [in] g The basic graph.
* @param [in] u The nodeset that induces the subgraph for which were checking the neighbourhood of v.
* @param [in] v A node, which might connect a circle in g[u].
* @returns True, if a neighbour of a neighbour of v is a neighbour of v.
*/
bool fvs::creates_circle(const Graph& g, const set<Node>& u, const Node& v) {
	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
	pair<edge_iterator, edge_iterator> eIt = out_edges(v, g);

	set<Node> neighbours;
	for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
		if (u.end() != u.find(target((*it), g))) {
			neighbours.insert(target((*it), g));
		}
	}

	for (const auto& i : neighbours) {
		eIt = out_edges(i, g);
		for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
			if (neighbours.end() != neighbours.find(target((*it), g))) {
				return true; //A neighbour of a neighbour of v is neighbour of v (in u).
			}
		}
	}

	return false;
}

/**
* @brief Finds a node in u, which has atleast two neighbours in v with respect to g.
*
* This function iterates over all elements of u to find a node which has two neighbours
* in v, where both u and v are sets of nodes of the same graph g. Returns a node, if one
* was found or a null_vertex() if no fitting node was found.
*
* @param [in] g The graph this is based on.
* @param [in] u A node set. 
* @param [in] v Another node set.
* @returns A node of u with atleast two neighbours in v, if such a node exists or a null_vertex()
*		otherwise.
*/
Node fvs::two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node> &v) {
	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;

	pair<edge_iterator, edge_iterator> eIt;
	int numNeighbours;

	for (const auto& i : u) {
		eIt = out_edges(i, g);
		for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
			if (v.end() != v.find(target((*it), g))) {
				++numNeighbours;
			}

			if (numNeighbours > 1) {
				return i;
			}
		}
	}
	return graph_traits<Graph>::null_vertex();
}

/**
* @brief Creates the induced subgraph g[u].
*/
void fvs::induced_subgraph(Graph &s, const Graph& g, const set<Node>& u) {
	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
	
	s.clear();
	for (const auto& i : u) {
		pair<edge_iterator, edge_iterator> eIt = out_edges(i, g);
		for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
			if (u.end() != u.find(target((*it), g))) {
				add_edge(i, target((*it), g), s);
			}
		}
	}
}

/**
* @brief Computes a feedback vertex set.
*
* This algorithm was proposed in the paper, that we have here:
*	https://www.dropbox.com/sh/ar26siyo2cjw6y1/AACqJWkA0YHXxkg5FTz2ZEeBa/ChenFLLV08_ImprovedAlgorithmsForFeedbackVertexSetProblems.pdf?dl=0
*
* This function right now resembles a direct copy of the pseudo code of the article.
*
* @param [in] g The graph to find a feedback vertex set for.
* @param [in] f A subgraph of g, for which f is a forest and g-f is a forest.
* @param [in] k The currently guessed size of a min. feedback vertex set.
* @returns A pair of a set of nodes and a bool. The set of nodes contains a part of the feedback
*		vertex set. The bool will be false, if the algorithm decides that there is no fvs.
*/
pair<set<Node>, bool> fvs::compute_fvs(Graph& g, set<Node> f, int k) {
	set<Node> fvs;
	pair<set<Node>, bool> retValue;

	if (k < 0 || (k == 0 && has_cycle(g))) {
		return make_pair(fvs, false);
	}

	if (!has_cycle(g)) {
		return make_pair(fvs, true);
	}

	Node w = two_neighbour_node(g, f); // A vertex of f which has least two neighbors in g-f.

	if (w != graph_traits<Graph>::null_vertex()) {
		if (creates_circle(g, g-f, w)) {
			Graph h(g);
			remove_vertex(w, h);
			f.erase(w);
			retValue = compute_fvs(h, f, k - 1);

			if (false == retValue.second) {
				return make_pair(fvs, false);
			}
			else {
				fvs = retValue.first;
				fvs.insert(w);
				return make_pair(fvs, true);
			}
		}
		else {
			Graph h(g);
			remove_vertex(w, h);
			f.erase(w);
			retValue = compute_fvs(h, f, k - 1);

			if (true == retValue.second) {
				fvs = retValue.first;
				fvs.insert(w);
				return make_pair(fvs, true);
			}
			else {
				f.erase(w);
				return compute_fvs(g, f, k);
			}
		}
	}
	else {
		w = get_lowest_degree_node(g, f);
		if (out_degree(w, g) < 2) {
			Graph h(g);
			remove_vertex(w, h);
			f.erase(w);
			retValue = compute_fvs(h, f, k - 1);
		} else {
			f.erase(w);
			return compute_fvs(g, f, k);
		}
	}
}

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
void fvs::read_graph(Graph&g, const char* filepath) {
	typedef graph_traits<Graph>::vertices_size_type NodeId;
	size_t invalid = (size_t)-1;

	ifstream file(filepath, ios::in);
	//Check if the file exists.
	if (!file.is_open())
	{
		throw std::runtime_error(("Cannot open file at %p.", filepath));
		return;
	}

	//Reset the graph.
	g.clear();

	NodeId num_nodes = 0;
	string line;

	//Get first line of the file.
	getline(file, line);

	//Convert to stringstream
	stringstream stream;
	stream << line;

	//Check if the stream is valid.
	if (!stream.good())
	{
		throw std::runtime_error(("Invalid file format in %p.", filepath));
	}

	//Control parameter.
	size_t num_edges = 0;

	//Edge parameter
	size_t source = invalid;
	size_t target = invalid;
	int weight = 0;

	//Initialize the graph.
	stream >> num_nodes;
	Graph g(num_nodes);

	//Set the control parameter.
	getline(file, line);
	stream.clear();
	stream << line;
	stream >> num_edges;

	stream.str("");
	stream.clear();

	//Read out all lines.
	while (getline(file, line))
	{
		stream.str(line);

		//Add the edge.
		stream >> source;
		stream >> target;
		add_edge(source, target, g);

		stream.clear();
	}

	//Check if the graph was read correctly.
	if (boost::num_edges(g) != num_edges)
	{
		throw std::runtime_error(("Error graph %p was not read correctly.", filepath));
	}
	return;
}
