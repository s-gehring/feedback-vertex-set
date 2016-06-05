
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
//	cout << "FINDING CYCLES" << endl;
  	bool b = false;
	CycleVisitor cv(b);
	depth_first_search<Graph, CycleVisitor>(g, visitor(cv));
//	cout << "Has cycle: " << (b ? "true" : "false") << endl;
	return b;
}


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
pair<set<Node>, bool> fvs::find_semidisjoint_cycle(const Graph& g)
{
	bool b = false;
	set<Node> sdc;
	SemiDisjointCycleVisitor cv(b ,sdc);
	depth_first_search<Graph, SemiDisjointCycleVisitor>(g, visitor(cv));
	pair<set<Node>, bool> retValue = make_pair(sdc, b);
	return retValue;
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
	typedef graph_traits<Graph>::vertices_size_type id;
	if (u.size() == 0) {
		return graph_traits<Graph>::null_vertex();
//		throw runtime_error("Error: Searching lowest degree node in an empty set.");
	}

	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
	pair<id, Node> lowestDegreeNode = make_pair(1000, graph_traits<Graph>::null_vertex());

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
* @param [in] u The nodeset that does not induce the subgraph for which were checking the neighbourhood of v.
* @param [in] v A node, which might connect a circle in g[u].
* @returns True, if a neighbour of a neighbour of v is a neighbour of v.
*/
bool fvs::creates_circle(const Graph& g, const set<Node>& u, const Node& v) {
	//Make sure we dont call the rest on an invalid call.
	if (v == graph_traits<Graph>::null_vertex()) {
		//Might be useful to find a proper way for this...
		return false; 
	}
	
/*	cout << "Checking for circles." << endl;
	cout << "Node that is checked: " << v << endl;
	cout << "Set: "; for (const auto& i : u) { cout << i << ", "; } cout << endl;*/

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
//				cout << "Circle found." << endl;
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
* @param [in] u A node subset of g which partitions g in u and g-u.
* @returns A node of u with atleast two neighbours in v, if such a node exists or a null_vertex()
*		otherwise.
*/
Node fvs::two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node> &v) {
	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
	typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;
	
/*	cout << "TWO NEIGHBOUR NODE:" << endl;
	cout << "Set 1: ";
	for (const auto& i : u) { cout << i << ", "; } cout << endl;
	cout << "Set 2: ";
	for (const auto& i : v) { cout << i << ", "; } cout << endl;
	print_graph(g);*/
	

	pair<edge_iterator, edge_iterator> eIt;

	for (const auto& i : u) {
		eIt = out_edges(i, g);
		
		set<Node> neighbours;
		pair<adjacency_iterator, adjacency_iterator> adj = adjacent_vertices(i, g);
		for (adjacency_iterator aIt = adj.first; aIt != adj.second; ++aIt) {
			neighbours.insert((*aIt));
		}
		
		for (const auto& n : neighbours) {
			pair<adjacency_iterator, adjacency_iterator> neighboursNeighbours = adjacent_vertices(n, g);
			for (adjacency_iterator aIt = neighboursNeighbours.first; aIt != neighboursNeighbours.second; ++aIt) {
				if (neighbours.find((*aIt)) != neighbours.end()) {
//					cout << "Node found: " << i << endl;
					return i;
				}
			}
		}
	}
//	cout << "No node found: " << endl;
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
* @brief Maintains integrity in the node set, s.t. the nodes remain correct references when another node is deleted.
*
* This method move all references of nodes in u, s.t. those references will still be correct, when `aDeletedNode` is
* deleted. The function assumes, that `u` does not contain `aDeletedNode`.
*
* References are updated in the same way that bgl uses when removing a vertex, all indices greater than `aDeletedNode`
* will be decremented by one.
*
* The direction_tag denotes wether this operation is performed upon deletion of a node (direction_tag::forward), or if
* the indices should be recomputed (when going up the recursion tree).
*
* @param [in] g A graph. (Might not be necessary).
* @param [in, out] u A set of nodes. It will be updated.
* @param [in] aDeletedNode The node that will be deleted soon.
* @param [in] dt A direction_tag, defaults to forward.
*/
void fvs::maintain_integrity(Graph& g, set<Node>& u, Node aDeletedNode, direction_tag dt) {
/*	cout << "Maintain integrity." << endl;
	cout << "Deleted Node: " << aDeletedNode << endl;
	cout << "U: ";				
	for (const auto& i : u) { cout << i << ", "; } cout << endl;*/
	
	set<Node> newSet;
	for (set<Node>::iterator it = u.begin(); it != u.end(); ++it) {
		if (direction_tag::forward == dt) { 
			if (*it > aDeletedNode && 0 < *it) {
				newSet.insert(((*it) - 1));
				//(*it) = (*it) - 1;
			}
			else if (*it < aDeletedNode || (0 == *it && 0 != aDeletedNode)) {
				newSet.insert(((*it)));
			}
		}
		else {
			if (*it >= aDeletedNode) {
				newSet.insert(((*it) + 1));
			}
			else if (*it < aDeletedNode) {
				newSet.insert((*it));
			}
		}
	}
	
	u = newSet;
	
	//cout << "New U: ";				
	//for (const auto& i : u) { cout << i << ", "; } cout << endl << endl;
}


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
pair<set<Node>, bool> fvs::compute_fvs(const Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k) {
	set<Node> fvs;
	pair<set<Node>, bool> retValue;
	
/*	cout << "Number of Nodes: " << num_vertices(g) << endl;
	cout << "k: " << k << endl;
	cout << "FVS: ";
	for (const auto& i : f) { cout << i << ", "; } cout << endl;
	cout << "V2: ";
	for (const auto& i : v2) { cout << i << ", "; } cout << endl;
	
	print_graph(g);*/

	if (k < 0 || (k == 0 && has_cycle(g))) {
//		cout << "Returning false, k got small and g has still cycles." << endl;
		return make_pair(fvs, false);
	}

	if (!has_cycle(g)) {
//		cout << "g has no cycle." << endl;
		//print_graph(g);
		return make_pair(fvs, true);
	}

	Node w = two_neighbour_node(g, f, v2); // A vertex of f which has least two neighbors in g-f.
//	cout << "w: " << w << endl;
	if (w != graph_traits<Graph>::null_vertex()) {
		if (creates_circle(g, v2, w)) {
//			cout << w << " creates a circle in g." << endl;
			
			Graph h(g);
			f.erase(w);
			clear_vertex(w, h);
			remove_vertex(w, h);
			maintain_integrity(h, f, w);
			maintain_integrity(h, v2, w);
			retValue = compute_fvs(orig, h, f, v2, k - 1);
			
//			cout << "Subcall retuned: " << (retValue.second ? "true" : "false") << endl;
//			cout << "Subcall-Set: "; for (const auto& i : retValue.first) { cout << i << ", "; } cout << endl;
			
			
			if (false == retValue.second) {
//				cout << "Returning false after subcall. " << endl;
				return make_pair(fvs, false);
			}
			else {
				fvs = retValue.first;
				
				maintain_integrity(h, fvs, w, direction_tag::inverse);
				
				fvs.insert(w);
				
//				cout << "Returning from 3.1: ";
//				for (const auto& i : fvs) { cout << i << ", "; } cout << endl;
				return make_pair(fvs, true);
			}
		}
		else {
			Graph h(g);
			f.erase(w);
			clear_vertex(w, h);
			remove_vertex(w, h);
			maintain_integrity(h, f, w);
			maintain_integrity(h, v2, w);
			retValue = compute_fvs(orig, h, f, v2, k - 1);

			if (true == retValue.second) {
				fvs = retValue.first;
				maintain_integrity(h, fvs, w, direction_tag::inverse);
				
				fvs.insert(w);
				
//				cout << "Returning: ";
//				for (const auto& i : f) { cout << i << ", "; } cout << endl;
				return make_pair(fvs, true);
			}
			else {
				v2.insert(w);
//				cout << "Removing " << w << " from the fvs." << endl;
				return compute_fvs(orig, g, f, v2, k);
			}
		}
	}
	else {
		w = get_lowest_degree_node(g, f);
		if (graph_traits<Graph>::null_vertex() != w && out_degree(w, orig) < 2) {
//			cout << "Out degree of " << w << " is smaller than 2" << endl;
			Graph h(g);
			f.erase(w);
			clear_vertex(w, h);
			remove_vertex(w, h);
			maintain_integrity(h, f, w);
			maintain_integrity(h, v2, w);

			return retValue = compute_fvs(orig, h, f, v2, k - 1);
		} 
		else if (w != graph_traits<Graph>::null_vertex()) {
			f.erase(w);
			v2.insert(w);

//			cout << "Removing " << w << " from the fvs." << endl;
			return compute_fvs(orig, g, f, v2, k);
		}
	}
	
	return make_pair(fvs, false);
}

pair<string, string> explode(string s) {
  return make_pair("", "");
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
	//size_t invalid = (size_t)-1;
  ifstream file(filepath, ios::in);
	//Check if the file exists.
	if (!file.is_open())
	{
		throw std::runtime_error("Cannot open file.");
		return;
	}

	//Reset the graph.
	g.clear();
	string line;
	//Read out all lines.
	//int lineNo = 0;
	

	while (getline(file, line))
	{
		istringstream iss(line);
	  int src, dst; 
	  int x = sscanf(line.c_str(), "%d %d", &src, &dst);
	  if(x != 2) {
	    throw std::runtime_error("Can't parse file. Wrong format?");
	    return;
	  }

	  // TODO: Add Lookup Table
		add_edge(vertex(src, g), vertex(dst, g), g);

	}
  file.close();
	return;
}

/**
* @brief: Deletes all vertices of degree at most 1 along with all incident edges from a given graph.
*
* As long as a given graph has vertices of degree at most 1, all incident edges and all vertices
* are deleted (by setting their weight to -1). We need this as a subroutine for the 2-approx-algo.
*
* @param [in] g The graph.
* @param [in] weights The weights of the vertices.
*/
void fvs::cleanup(Graph& g, map<Node, double>& weights)
{
	graph_traits<Graph>::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end && num_vertices(g) > 0; ++vi)
	{
		if (in_degree(*vi, g) <= 1 && weights[*vi] >= 0)
		{
			// remove edges and delete
			clear_vertex(*vi, g);
			weights[*vi] = -1;
			// check everything again
			tie(vi, vi_end) = vertices(g);
			// we can definitively save some runtime by just checking the neighbours we have already visited
		}
	}
}

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
set<Node> fvs::two_approx_fvs(const Graph& orig)
{
	// initialize weights and f
	Graph g(orig);
	map<Node, double> weights;
	graph_traits<Graph>::vertex_iterator vi, vi_end;
	set<Node> f;
	stack<Node> s;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		weights[*vi] = in_degree(*vi, g) - 1;
		// maybe there is a smarter way to initialize the weights, this seems to be reasonable fast due to its degree-proportionality
	}
	cleanup(g, weights);
	// count number of vertices with weight > 0
	int active_vertices = 0;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
		if (weights[*vi] >= 0) {
			active_vertices++;
		}
	}
	while (active_vertices > 0)
	{
		// contains a semidisjoint cycle?
		// Here, we do not the same as the algorithm in the paper: we only put the node with degree >2 into the fvs
		// and delete all the others. In O notation this is no difference but in practice we safe a tiny bit of time.
		pair<set<Node>, bool> sdcycle = find_semidisjoint_cycle(g);
		if (sdcycle.second) {
			// find max element
			bool true_cycle = true;
			for (set<Node>::iterator it = sdcycle.first.begin(); it != sdcycle.first.end(); ++it) {
				if (in_degree(*it, g) > 2) {
					weights[*it] = 0; // add to f
					true_cycle = false;
					sdcycle.first.erase(it);
					it = sdcycle.first.end(); // stop searching
					it--;
				}
			}
			if (true_cycle) {
				weights[*sdcycle.first.begin()] = 0; // just any of them will do the job
			}
			// delete all the others
			for (set<Node>::iterator it = sdcycle.first.begin(); it != sdcycle.first.end(); ++it) {
				if (it == sdcycle.first.begin() && !true_cycle) {
					clear_vertex(*sdcycle.first.begin(), g);
					weights[*sdcycle.first.begin()] = -1;
				}
				else {
					clear_vertex(*it, g);
					weights[*it] = -1;
				}
			}
		}
		else { // is clean and contains no semidisjoint cycle
			double gamma = 0;
			for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
			{
				if ((weights[*vi] / (in_degree(*vi, g) - 1) < gamma || gamma == 0) && weights[*vi] >=0) {
					gamma = weights[*vi]/(in_degree(*vi, g)-1);
				}
			}
			for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
				if (weights[*vi] >= 0) {
					weights[*vi] = weights[*vi] - gamma*(in_degree(*vi, g)-1);
				}
			}
		}
		// handle vertices with weight 0
		for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
			if (weights[*vi] == 0) {
				f.insert(*vi);
				s.push(*vi);
				clear_vertex(*vi, g);
				weights[*vi] = -1; // this means that the vertex is considered to be deleted
			}
		}
		cleanup(g, weights);
		// count number of vertices with weight > 0
		active_vertices = 0;
		for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
			if (weights[*vi] >= 0) {
				active_vertices++;
			}
		}
	}
	// shrink the approximation of the fvs set
	while (!s.empty()) {
		Node u = s.top();
		s.pop();
		// construct graph to be checked
		Graph g(orig);
		set<Node>::iterator it_u;
		for (set<Node>::iterator it = f.begin(); it != f.end(); ++it) {
			if (*it != u) {
				clear_vertex(*it, g); // the vertex does not need to be deleted
			}
			else {
				it_u = it;
			}
		}
		// is f without u an fvs in the original g?
		if (!has_cycle(g)) {
			f.erase(it_u);
		}
	}
	return f;
}

/**
* @brief Print g to stdout.
*/
void fvs::print_graph(const Graph& g) {
	typedef graph_traits<Graph>::vertex_iterator node_iterator;
	typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
	
	cout << "Printing a graph." << endl;
	cout << "Number of nodes: " << num_vertices(g) << endl;
	cout << "Number of edges: " << num_edges(g) << endl;
	if(num_edges(g) > 1000 || num_vertices(g) > 500) return;
	pair<node_iterator, node_iterator> nIt = vertices(g);
	for (node_iterator it = nIt.first; it != nIt.second; ++it) {
		cout << "Edges outgoing from " << (*it) << ":" << endl;
		pair<edge_iterator, edge_iterator> eIt = out_edges((*it), g);
		for (edge_iterator edgeIt = eIt.first; edgeIt != eIt.second; ++edgeIt) {
			cout << source((*edgeIt), g) << " -> " << target((*edgeIt), g) << endl;
		}
	}
	cout << "---------------------------" << endl;
}
