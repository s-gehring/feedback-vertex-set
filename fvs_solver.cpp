
#include <utility>
#include <limits.h>
#include <fstream>
#include <sstream>

#include "fvs_solver.h"

using namespace fvs;

/**
* @brief Checks, wether the given graph contains a cycle.
*
* This function uses dfs and some shortcuts to check if the given graph contains cycles.
*
* @param [in] g The graph.
* @returns True, if there is a cycle in g.
*/
bool fvs::has_cycle(Graph& g) {
	return g.has_cycle();
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
pair<set<Node>, bool> fvs::find_semidisjoint_cycle(Graph& g)
{
	return g.find_semidisjoint_cycle();
}

/**
* This is basically an example on how to do that with boost.
* Dont use this method, just call the original.
*/
bool fvs::edge_exists_between(Graph& g, Node u, Node v) {
	return g.has_edge(u,v);
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
Node fvs::get_lowest_degree_node(Graph &g, const set<Node>& u) {
	return g.lowest_deg_node(u);
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
bool fvs::creates_circle(Graph& g, const set<Node>& u, const Node& v) {
    if(v == INVALID_NODE) {
        cout << "Warning: Called creates_circle() with invalid node." << endl;
        return false;
    }
/*	cout << "Checking for circles." << endl;
	cout << "Node that is checked: " << v << endl;
	cout << "Set: "; for (const auto& i : u) { cout << i << ", "; } cout << endl;*/

	pair<Neighborhood, bool> neighbors = g.get_neighbors(v);
    if(!neighbors.second) {
        cout << "Warning: Called creates_circle() with a vertex (v) not being in the graph (g)."<<endl;
        return false;
    }
    
	set<Node> neighbors_in_u;
	for(const auto &it : neighbors.first) {
	    if(u.find(it) != u.end()) {
	        neighbors_in_u.insert(it);
	    }
	}

	for (const auto& i : neighbors_in_u) {
		auto eIt = g.get_neighbors(i); // Neighborhood , bool
		if(!eIt.second) {
		    throw std::runtime_error("Error: There is a vertex which is in u, but not in g.");
		    return false;
		}
		for (const auto& it : eIt.first) {
		    if(neighbors_in_u.find(it) != neighbors_in_u.end()) {
		        return true;
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
Node fvs::two_neighbour_node(Graph& g, const set<Node> &u, const set<Node> &v) {
	
	for (const auto& i : u) {
		int neighbors = 0;
		Neighborhood n = g.get_neighbors(i).first;
		for(const auto& j : n) {
		    if(v.find(j) != v.end()) {
		        if(++neighbors > 1) return i;
		    }
		}
	}
	return INVALID_NODE;
}

/**
* @brief Creates the induced subgraph g[u].
*/
void fvs::induced_subgraph(Graph &s, Graph& g, const set<Node>& u) {
	s.clear();
	for (const auto& i : u) {
		std::pair<Neighborhood, bool> eIt = g.get_neighbors(i);
		if(!eIt.second) {
		    cout << "Warning: Called induced_subgraph with a nodeset containing at least one node not in g"<<endl;
		} else {
		    for(const auto& j : eIt.first) {
		        if(u.find(j)!=u.end()) {
		            s.add_edge(i, j);
		        }
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
    return; // Obsolete
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
pair<set<Node>, bool> fvs::compute_fvs(Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k) {
	set<Node> fvs;
	pair<set<Node>, bool> retValue;
	
	if (k < 0 || (k == 0 && has_cycle(g))) {
		//cout << "Returning false, k got small and g has still cycles." << endl;
		return make_pair(fvs, false);
	}

	if (!has_cycle(g)) {
	//	cout << "g has no cycle." << endl;
		//print_graph(g);
		return make_pair(fvs, true);
	}

	Node w = two_neighbour_node(g, f, v2); // A vertex of f which has least two neighbors in g-f.
	//cout << "w: " << w << endl;
	if (w != INVALID_NODE) {
		if (creates_circle(g, v2, w)) {
			//cout << w << " creates a circle in g." << endl;
			
			Graph h(g);
			f.erase(w);
			h.remove_node(w);
			retValue = compute_fvs(orig, h, f, v2, k - 1);
			
		//	cout << "Subcall retuned: " << (retValue.second ? "true" : "false") << endl;
		//	cout << "Subcall-Set: "; for (const auto& i : retValue.first) { cout << i << ", "; } cout << endl;
			
			
			if (false == retValue.second) {
		//		cout << "Returning false after subcall. " << endl;
				return make_pair(fvs, false);
			}
			else {
				fvs = retValue.first;
				fvs.insert(w);
				
		//		cout << "Returning from 3.1: ";
				//for (const auto& i : fvs) { cout << i << ", "; } cout << endl;
				return make_pair(fvs, true);
			}
		}
		else {
			Graph h(g);
			f.erase(w);
			h.remove_node(w);
			retValue = compute_fvs(orig, h, f, v2, k - 1);

			if (true == retValue.second) {
				fvs = retValue.first;
				
				fvs.insert(w);
				
//cout << "Returning: ";
			//	for (const auto& i : f) { cout << i << ", "; } cout << endl;
				return make_pair(fvs, true);
			}
			else {
				v2.insert(w);
	//			cout << "Removing " << w << " from the fvs." << endl;
				return compute_fvs(orig, g, f, v2, k);
			}
		}
	}
	else {
		w = get_lowest_degree_node(g, f);
		if (INVALID_NODE != w && orig.get_single_degree(w) < 2) {
//			cout << "Out degree of " << w << " is smaller than 2" << endl;
			Graph h(g);
			f.erase(w);
      h.remove_node(w);

			return retValue = compute_fvs(orig, h, f, v2, k - 1);
		} 
		else if (w != INVALID_NODE) {
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
void fvs::read_graph(Graph &g, const char* filepath) {
    ifstream file(filepath, ios::in);
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
		g.add_edge(src, dst);

	}
    file.close();
	return;
}

/**
* @brief: Deletes all vertices of degree at most 1 along with all incident edges from a given graph.
*
* As long as a given graph has vertices of degree at most 1, all incident edges and all vertices
* are deleted. We need this as a subroutine for the 2-approx-algo.
*
* @param [in] g The graph.
*/
void fvs::cleanup(Graph& g)
{
	Node neighbour;
	Node help;
	for(const auto &it : g.get_adjacency_list()) {
		if(it.second.size() == 0) {
			g.remove_node(it.first);
		}
		else if (it.second.size() == 1) {
			neighbour = *(it.second.begin());
			g.remove_node(it.first);
			// check all the neighbouring nodes we have already checked
			while (g.get_single_degree(neighbour) < 2 && neighbour < it.first) {
				if (g.get_single_degree(neighbour) == 0) {
					g.remove_node(neighbour);
					neighbour = it.first; // stop checking neighbours
				}
				else if (g.get_single_degree(neighbour) == 1) {
					help = neighbour;
					neighbour = *(g.get_neighbors(help).first.begin()); // for repeating the process
					g.remove_node(help);
				}
			}
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

set<Node> fvs::two_approx_fvs(Graph& orig)
{
	Graph g(orig); // our working copy
	set<Node> f; // the approx of the fvs
	map<Node, double> weights; // individual weights for each node
	stack<Node> s; // here we save all nodes of the fvs for checking their neccessarity in the end
	cleanup(g);
	// initialize weights
	for (const auto &it : g.get_adjacency_list()) {
		weights[it.first] = g.get_single_degree(it.first) - 1;
		// maybe there is a faster way to do it but this is reasonable fast due to its degree proportionality
	}
	while (g.n > 0) {
		// contains a semidisjoint cycle?
		// Here, we do not the same as the algorithm in the paper: we only put the node with degree >2 into the fvs
		// and delete all the others. In O notation this is no difference but in practice we safe a tiny bit of time.
		pair<set<Node>, bool> sdcycle = find_semidisjoint_cycle(g);
		if (sdcycle.second) {
			// find max element
			bool true_cycle = true;
			for (set<Node>::iterator it = sdcycle.first.begin(); it != sdcycle.first.end(); ++it) {
				if (g.get_single_degree(*it) > 2) {
				
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
					g.remove_node(*sdcycle.first.begin());
				}
				else {
					g.remove_node(*it);
				}
			}
		} else { // is clean and contains no semidisjoint cycle
			double gamma = 0;
			// find minimum
			for (const auto &it : g.get_adjacency_list()) {
				if ((weights[it.first] / (g.get_single_degree(it.first)-1) < gamma || gamma == 0)) {
					gamma = weights[it.first]/(g.get_single_degree(it.first)-1);
				}
			}
			// update weights
			for (const auto &it : g.get_adjacency_list()) {
				weights[it.first] = weights[it.first] - gamma*(g.get_single_degree(it.first)-1);
			}
		}
		// handle remaining vertices with weight 0
		for (const auto &it : g.get_adjacency_list()) {
			if (weights[it.first] == 0) {
				f.insert(it.first);
				s.push(it.first);
				g.remove_node(it.first);
			}
		}
		cleanup(g);
	}
	// shrink the approximation of the fvs set
	//Graph h(orig);
	stack<Node> nodes_to_reset;
	stack<Edge> edges_to_reset;
	Node u;
	set<Node>::iterator it_u;
	Neighborhood neighbours;
	while (!s.empty()) {
		Graph h(orig); // # copy whole graph
		u = s.top();
		s.pop();
		// construct graph to be checked
		for (set<Node>::iterator it = f.begin(); it != f.end(); ++it) {
			if (*it != u) {
				/*
				// save Node and all its incident edges, to be able to reset it again
				nodes_to_reset.push(*it);
				neighbours = h.get_neighbors(*it).first;
				for (Neighborhood::iterator it1 = neighbours.begin(); it1 != neighbours.end(); ++it1) {
					edges_to_reset.push(make_pair(*it, *it1));
				}
				*/
				h.remove_node(*it);
			}
			else {
				it_u = it;
			}
		}
		// is f without u an fvs in the original g?
		if (!has_cycle(h)) {
			f.erase(it_u);
		}
		/*
		// reset graph again
		while (!nodes_to_reset.empty()) {
			h.add_node(nodes_to_reset.top());
			nodes_to_reset.pop();
		}
		while (!edges_to_reset.empty()) {
			h.add_edge(edges_to_reset.top().first, edges_to_reset.top().second);
			edges_to_reset.pop();
		}
		*/
	}
	return f;
}

/**
* @brief Print g to stdout.
*/
void fvs::print_graph(Graph& g) {
	
	cout << "Printing a graph." << endl;
	cout << "Number of nodes: " << g.n << endl;
	cout << "Number of edges: " << g.m << endl;
	if(g.m > 1000 || g.n > 500) {
	  cout << "Graph too big, skipping complete printing."<<endl;
	  return;
	}
	for (const auto &it : g.get_adjacency_list()) {
		cout << "Edges outgoing from " << it.first << ":" << endl;
		for (const auto &eit : it.second) {
			cout << it.first << " -> " << eit << endl;
		}
	}
	cout << "---------------------------" << endl;
}
