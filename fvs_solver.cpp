#include <utility>
#include <limits.h>
#include <fstream>
#include <sstream>

#include "fvs_solver.hpp"

using namespace fvs;

bool fvs::has_cycle(Graph& g) {
	return g.has_cycle();
}

pair<list<Node>, bool> fvs::find_semidisjoint_cycle(Graph& g)
{
	return g.find_semidisjoint_cycle();
}

Node fvs::get_lowest_degree_node(Graph &g, const set<Node>& u) {
	return g.lowest_deg_node(u);
}

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

void fvs::cleanup(Graph& g)
{
	Node neighbour;
	Node help;
	set<Node> processed;
	for(auto &it : g.get_low_degree_nodes()) {
		if(g.get_single_degree(it) == 0) {
			g.remove_node(it);
		}
		else if (g.get_single_degree(it) == 1) {
			// check the neighbour and his neighbours if they were already processed
			neighbour = *(g.get_neighbors(it).first.begin());
			g.remove_node(it);
			while (g.get_single_degree(neighbour) < 2 && processed.find(neighbour) != processed.end()) {
				if (g.get_single_degree(neighbour) == 0) {
					g.remove_node(neighbour);
					neighbour = it; // stop checking neighbours
				}
				else if (g.get_single_degree(neighbour) == 1) {
					help = neighbour;
					neighbour = *(g.get_neighbors(help).first.begin()); // for repeating the process
					g.remove_node(help);
				}
			}
		}
		processed.insert(it);
	}
}

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
				g.remove_node(*it);
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

bool fvs::is_fvs(const Graph& g, const set<Node>& fvs)
{
	Graph h(g);
	for (set<Node>::iterator it = fvs.begin(); it != fvs.end(); ++it) {
		h.remove_node(*it);
	}
	return !has_cycle(h);
}

pair<set<Node>, bool> fvs::compression_fvs(const Graph& orig, const set<Node>& S) {
	Graph g(orig);
	int k = S.size() - 1;
	int n = pow(2, k);
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
			set<Node> v_without_s = set_minus(V, S);
			set<Node> s_without_d = set_minus(S, D);
			result = forest_bipartition_fvs(g,g,v_without_s,s_without_d,k-D.size());
			if (result.second) {
				return make_pair(set_union(result.first, D), true);
			}
		}
		D.clear();
	}
	return make_pair(S, false);
}

set<Node> fvs::set_minus(const set<Node> S, const set<Node> T) {
	set<Node> difference;
	for (set<Node>::iterator it = S.begin(); it != S.end(); ++it){
		if (T.find(*it) == T.end()) {
			difference.insert(*it);
		}
	}
	return difference;
}

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
