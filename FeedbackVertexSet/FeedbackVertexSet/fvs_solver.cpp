
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
* PLEASE SOMEBODY COME UP WITH A NAME FOR U. A DESCRIPTIVE ONE!
*
* @param g [in] The graph.
* @param u [in] A set of nodes. 
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