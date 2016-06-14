
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
pair<list<Node>, bool> fvs::find_semidisjoint_cycle(Graph& g)
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
pair<set<Node>, bool> fvs::forest_bipartition_fvs(Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k) {
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
			retValue = forest_bipartition_fvs(orig, h, f, v2, k - 1);
			
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
			retValue = forest_bipartition_fvs(orig, h, f, v2, k - 1);

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
				return forest_bipartition_fvs(orig, g, f, v2, k);
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

			return retValue = forest_bipartition_fvs(orig, h, f, v2, k - 1);
		} 
		else if (w != INVALID_NODE) {
			f.erase(w);
			v2.insert(w);

//			cout << "Removing " << w << " from the fvs." << endl;
			return forest_bipartition_fvs(orig, g, f, v2, k);
		}
	}
	
	return make_pair(fvs, false);
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
	set<Node> processed;
	for(auto &it : g.get_adjacency_list()) {
		if(g.get_single_degree(it.first) == 0) {
			g.remove_node(it.first);
		}
		else if (g.get_single_degree(it.first) == 1) {
			// check the neighbour and his neighbours if they were already processed
			neighbour = *(g.get_neighbors(it.first).first.begin());
			g.remove_node(it.first);
			while (g.get_single_degree(neighbour) < 2 && processed.find(neighbour) != processed.end()) {
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
		processed.insert(it.first);
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
		weights[it.first] = 1; // weights all need to be the same s.t. nodes with high degree are preferred later
	}
	while (g.get_n() > 0) {
		// contains a semidisjoint cycle?
		// Here, we do not the same as the algorithm in the paper: we only put the node with degree >2 into the fvs
		pair<list<Node>, bool> sdcycle = find_semidisjoint_cycle(g);
		if (sdcycle.second) {
			// take high degree element of semidisjoint cycle, if its a true disjoint cycle any element is fine, thus this will always do the job
			weights[sdcycle.first.front()] = 0;
			sdcycle.first.pop_front();
			// delete all other elements from the graph
			for (list<Node>::const_iterator it = sdcycle.first.begin(); it != sdcycle.first.end(); ++it)
			{
			    /**
			    ** Don't change the set you're iterating over!
			    ** You can invoke sdcycle.clear() after this loop
			    ** terminates, but keep your
          }
   filthy hands off
			    ** the list while you're using its iterators!
			    */
				g.remove_node(*it);
				// <evil>
    			//	sdcycle.first.erase(it);
    			// </evil>
			
			}
		}
		else { // is clean and contains no semidisjoint cycle
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
	Node u;
	while (!s.empty()) {
		u = s.top();
		f.erase(u);
		s.pop();
		// insert u again into the fvs if we have cycles without it
		if (!is_fvs(orig,f)) {
			f.insert(u);
		}
	}
	return f;
}

/**
*@brief: Checks whether a given set is a feedback vertex set in a given graph.
*
* @param[in] g The graph.
* @param[in] fvs The set to be checked.
* @returns True, if the given set is a feedback vertex set in g.
*/
bool fvs::is_fvs(const Graph& g, const set<Node>& fvs)
{
	// copy graph
	Graph h(g);
	for (set<Node>::iterator it = fvs.begin(); it != fvs.end(); ++it) {
		h.remove_node(*it);
	}
	return !has_cycle(h);
}

/**
*@brief: Tries to decrease the size of a given fvs in a given graph by 1.
*
* @param[in] orig The graph.
* @param[in] S The feedback vertex set which will be compressed.
* @returns A new fvs with size of the old fvs -1 if it exists, otherwise false.
*/
pair<set<Node>, bool> fvs::compression_fvs(const Graph& orig, const set<Node>& S) {
	Graph g(orig);
	int k = S.size() - 1;
	int n = pow(2, S.size() - 1);
	set<Node> D; // the guessed intersection
	// get nodes of the graph
	set<Node> V;
	for (const auto &it : g.get_adjacency_list()) {
		V.insert(it.first);
	}
	bool current;
	// convert set S to a vector
	vector<Node> T;
	for (set<Node>::iterator it = S.begin(); it != S.end(); ++it) {
		T.push_back(*it);
	}
	pair<set<Node>, bool> result;
	int h;
	for (int j = 0; j < n; j++) {
		// guess intersection using binary coding
		h = j;
		for (int l = 0; l < k; l++) {
			current = h % 2;
			if (current) {
				D.insert(T[l]);
			}
			h = (h - current) / 2;
		}
		set<Node> H = set_minus(S, D);
		Graph g(orig);
		for (const auto &it : g.get_adjacency_list()) {
			if (H.find(it.first) == H.end()) {
				g.remove_node(it.first);
			}
		}
		if (!has_cycle(g)) {
			// compute G without D
			Graph g(orig);
			for (set<Node>::iterator it = D.begin(); it != D.end(); ++it) {
				g.remove_node(*it);
			}
			// Invalid Initialization of non-const reference of type set.
			// This is why we have to save the two set_minus :/
			// However, I just called these sets v1 and v2, this is
			// an arbritrary choice I made and these names
			// don't correspond to any meaning of any paper.
			set<Node> v1 = set_minus(V, S);
			set<Node> v2 = set_minus(S, D);
			
			result = forest_bipartition_fvs(g,g,v1,v2,k-D.size());
			if (result.second) {
				return make_pair(set_union(result.first, D), true);
			}
		}
		D.clear();
	}
	return make_pair(S, false);
}

/**
*@brief: Computes the difference set of two given sets.
*
* @param[in] S The first set.
* @param[in] T The set the will be substracted from the first one.
* @returns The difference of the sets.
*/
set<Node> fvs::set_minus(const set<Node> S, const set<Node> T) {
	set<Node> difference;
	for (set<Node>::iterator it = S.begin(); it != S.end(); ++it){
		if (T.find(*it) == T.end()) {
			difference.insert(*it);
		}
	}
	return difference;
}

/**
*@brief: Computes the union of two given sets.
*
* @param[in] S The first set.
* @param[in] T The second set.
* @returns The union of the sets.
*/
set<Node> fvs::set_union(const set<Node> S, const set<Node> T) {
	set<Node> uni;
	set<Node> smaller_set;
	if (S.size() <= T.size()) {
		uni = T;
		smaller_set = S;
	}
	else {
		uni = S;
		smaller_set = T;
	}
	for (set<Node>::iterator it = smaller_set.begin(); it != smaller_set.end(); ++it) {
		uni.insert(*it);
	}
	return uni;
}

/**
* @brief Print g to stdout.
*/
void fvs::print_graph(Graph& g) {
	
	cout << "Printing a graph." << endl;
	cout << "Number of nodes: " << g.get_n() << endl;
	cout << "Number of edges: " << g.get_m() << endl;
	if(g.get_m() > 1000 || g.get_n() > 500) {
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
