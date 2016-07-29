
#include "fvs_solver.hpp"

using namespace fvs;
using namespace FvsGraph;
using namespace BinCount;
#ifndef debug
    #define debug if(0)
#endif
#define MAX_OUTPUT_ARTICULATION 40

set<Node> greedy(Graph g)
{
	set<Node> result;
	Graph h(g);
	h.delete_low_degree_nodes();
	while(h.get_n()!=0)
	{
		Node n=-1;
		int maxDegree=-1;
    	for (const auto &it : h.get_adjacency_list()) {
			int degree=h.get_neighbors(it.first).first.size();
			if(maxDegree<degree)
			{
				maxDegree=degree;
				n=it.first;
			}
		}
		result.insert(n);
		h.remove_node(n);
		h.delete_low_degree_nodes();
	}
	return result;
}

/*
*  Output all articulation vertices and bridges.
*  Note, that I defined that
*  "An element is an articulation element iff
*  it is an articulation node or a bridge."
*
*  Basically something we can cut the graph at
*  to simplify computation (hopefully by
*  an exponential factor).
*
* @param [in] g The graph used.
*/
void output_arti_elems_bridges(const Graph& g) {
	pair<set<Node>, unordered_set<Edge> > master_of_arts = g.get_articulation_elements();
	/*
	**  Output all artiulcation elements.
	*/
	cout << "Articulation Vertices:" << endl << "{";
	bool first_out = true;
	int max_output = MAX_OUTPUT_ARTICULATION;
	for (const auto &it : master_of_arts.first) {
		if (--max_output == 0) { cout << "...[" << master_of_arts.first.size() - MAX_OUTPUT_ARTICULATION << " more]"; break; }
		if (first_out) {
			first_out = false;
			cout << it;
		}
		else {
			cout << ", " << it;
		}
	}
	cout << "}" << endl << endl;
	cout << "Bridges:" << endl << "{";
	first_out = true;
	max_output = MAX_OUTPUT_ARTICULATION;
	for (const auto &it : master_of_arts.second) {
		if (--max_output == 0) { cout << "...[" << master_of_arts.second.size() - MAX_OUTPUT_ARTICULATION << " more]"; break; }
		if (first_out) {
			first_out = false;
			cout << "(" << it.first << "," << it.second << ")";
		}
		else {
			cout << ", (" << it.first << "," << it.second << ")";
		}
	}
	cout << "}" << endl << endl;
}

/*
* Runs iterative compression to find a minimum feedback vertex set in a given graph.
*
* @param [in] g The filepath of the file where the graph is in.
*/
set<Node> run_iter_comp(Graph g) {
    // compute min fvs by using iterative compression
    set<Node> min_fvs = compute_min_fvs(g);
    // sanity check and output of the results
    if (is_fvs(g, min_fvs)) {
	return min_fvs;
    }
    else {
    debug cout << "Error: The set we have found is not an FVS!" << endl;
    return set<Node>();
    }
}

void remove_bridges(Graph &g, const unordered_set<Edge> &bridges) {
    //debug  cout << "Removing edges: ";
    for(const auto &it : bridges) {
      //debug cout << "(" << g.get_node_name(it.first) << "," << g.get_node_name(it.second)<< "),";
      g.remove_edge(it.first, it.second);
      if(g.get_single_degree(it.first) == 0) {
          g.remove_node(it.first);
      }
      if(g.get_single_degree(it.second) == 0) {
          g.remove_node(it.second);
      }
    }
    //debug cout <<endl<<endl;
    debug cout << "Deleted " << bridges.size() << " bridges (useless edges)." <<endl;
}

set<Edge> contract_edges(Graph &g) {
    set<Edge> branching_pairs;
    Node u = INVALID_NODE;
    Node v = INVALID_NODE;
    int contracted_edges = 0;
    //debug cout << "Contracting edges: ";
    unordered_set<Node> to_delete;
    unordered_set<Node> candidates = g.get_low_degree_nodes();
    for(auto &it : candidates) {
        if(g.get_single_degree(it) == 2) {
            // Shitty hack, because iterators are weird:
            bool first = true;
            for(const auto &it2 : g.get_neighbors(it).first) {
                if(first) {
                    first = false;
                    u = it2;
                } else {
                    v = it2;
                }
            }
	    /*
	    * If we want to contract this edge, we would create multiedges
	    * instead, we keep a single edge and remember that we have to take one of them into the fvs
	    */
            if(g.has_edge(u, v)) {
		branching_pairs.insert(make_pair(u, v));
	        to_delete.insert(it); // do not delete it now, to make sure that the multiedge nodes are not deleted
            }
	    else {
		g.add_edge(u, v);
	        g.remove_node(it);
	    }
	    ++contracted_edges;
        }
    }
    for (auto &it : to_delete) {
	g.remove_node(it);
    }
    debug cout << "Contracted " << contracted_edges << " edges." <<endl;
    return branching_pairs;
}

/*
* We basically try to find all minimal vertex covers in the graph only consisting of the multiedges.
* Since #multiedges < k it is okay to do this recursive branching.
*/
set<set<Node>> multi_edge_partitions(set<set<Node>>& m, set<Node>& taken,  Graph& g) {
	set<set<Node>> return_value;
	AdjacencyList adj = g.get_adjacency_list();
	// graph is partitioned
	if (adj.size() == 0) {
		return_value = m;
		return_value.insert(taken);
		return return_value;
	}
	// branch on one node
	else {
		Node branched_node = adj.begin()->first;
		// 1) take it
		Graph h1(g);
		set<Node> new_taken = taken;
		new_taken.insert(branched_node);
		// do not take neighbours
		// only remove the edges otherwise other multiedges are destroyed too
		// this is simply done by removing the node itself
		h1.remove_node(branched_node);
		// but we can delete neighbours that are now left alone
		for (auto &it : adj.begin()->second) {
			if (h1.get_single_degree(it) == 0) {
				h1.remove_node(it);
			}
			// or if it has degree 1 and again a neighbour with degree 0
			// we can skip it and take its neighbour
			else if (h1.get_single_degree(it) == 1) {
				if (h1.get_single_degree(*h1.get_neighbors(it).first.begin()) == 0) {
					new_taken.insert(*h1.get_neighbors(it).first.begin());
					h1.remove_node(*h1.get_neighbors(it).first.begin());
					h1.remove_node(it);
				}
			}
		}
		set<set<Node>> new_partition1 = multi_edge_partitions(m, new_taken, h1);
		m.insert(new_partition1.begin(), new_partition1.end());
		// 2) do not take it
		Graph h2(g);
		new_taken = taken;
		// take all neighbours
		for (auto &it : adj.begin()->second) {
			new_taken.insert(it);
			h2.remove_node(it);
		}
		// and delete all neighbouring nodes that now have degree 0
		for (auto &it : adj) {
			if (h2.get_single_degree(it.first) == 0 && it.first != branched_node) {
				h2.remove_node(it.first);
			}
		}
		h2.remove_node(branched_node);
		set<set<Node>> new_partition2 = multi_edge_partitions(m, new_taken, h2);
		m.insert(new_partition2.begin(), new_partition2.end());
	}
	return m;
}

int main(int argc, char** argv) {
	long seed = 0;
    if(argc > 1) {
        // Guess that the last string is the seed.
        // For everything else, we can define undefined behavior
        // and just not care.
        char* sd = argv[argc-1];
        if(sscanf(sd, "%li", &seed) == 0) {
            seed = 0;
            debug cout << "Can't parse seed. Setting to 0."<<endl;
        } else {
            debug cout << "Settings seed to "<<seed<<"."<<endl;   
        }
    }
    
	Galois & ga = Galois::getInstance();
	ga.set_w(16);
	ga.set_mode_logtb();
    ga.seed(seed);
	
	GraphData graph_data = read_graph();
    /*
    **  Edge case: n = 0
    */
    if(graph_data.graph.get_n() == 0) {
       debug cout << "Edge case: No not-necessary vertices. Outputting necessary vertices only."<<endl;
       for(const auto &it : graph_data.necessary_nodes) {
          cout << graph_data.mapping.second[it] << endl;
       }
       return 0;
    }
    
    
    
    
	set<Node> necessary_nodes = graph_data.necessary_nodes;
	Mapping node_names = graph_data.mapping;
	Graph g = graph_data.graph;

	g.assign_names(node_names);
	
	// Create original copy.
	Graph orig(g);

	// preprocessing
	remove_bridges(g, g.get_articulation_elements().second);

	//iteratively remove semidisjoint cycles
	bool progress = true;
	int sd_vertices = 0;
	while (progress) {
		pair<list<Node>, bool> sdcycle = find_semidisjoint_cycle(g);
		if (sdcycle.second) {
			necessary_nodes.insert(sdcycle.first.front());
			g.remove_node(sdcycle.first.front());
			sdcycle.first.pop_front();
			g.delete_low_degree_nodes();
			sd_vertices++;
		}
		progress = sdcycle.second;
	}
	debug cout << "Found " << sd_vertices << " vertices of the FVS while removing semidisjoint cycles." << endl;
	
	// Split graph on connected components.  
	std::list<Graph> connected_graphs;
	get_connected_graphs(g, g.get_connected_components(), connected_graphs);
	list<set<Node> > partial_solutions;
	set<Node> complete_solution;
	for (auto &it : connected_graphs) {
	  debug cout << "Trying to find partial solution for graph " << it.get_name() << "/"<<connected_graphs.size()<<" [n=" << it.get_n() << "|m=" << it.get_m() << "]" << endl;

	  set<Node> partial_solution;
	  if(degree3 && it.is_deg_three()) {
      debug cout << "Degree 3 case" << endl;
      set<Node> V1;
      int maxIndex=-1;
      for (const auto &it2 : it.get_adjacency_list())
      {
          V1.insert(it2.first);
          if (it2.first>maxIndex)
          {
              maxIndex=it2.first;
          }
      }
      vector<int> nodeToComponent(maxIndex+1);
      partial_solution = solveDegree3(it,V1,nodeToComponent);
		}
		else {
			// contract edges and start branching on multiedges within one connected component
			set<Edge> branching_pairs = contract_edges(it);
			if (branching_pairs.size() > 0) {
				debug cout << "Start branching on " << branching_pairs.size() << " different multi-edges ";
				debug cout << "for one connected component." << endl;
				// create graph by using all multi-edges
				Graph m;
				m.add_edges(branching_pairs);
				// generate all partitions
				set<Node> taken;
				set<set<Node>> partitions;
				partitions = multi_edge_partitions(partitions, taken, m);
				debug cout << "Reduced number of necessary branchings on multiedges from " << pow(2, branching_pairs.size());
				debug cout << " different branchings to " << partitions.size() << "." << endl;
				for (set<set<Node>>::iterator partition = partitions.begin(); partition != partitions.end(); ++partition) {
					Graph h(it);
					set<Node> current_solution;
					for (set<Node>::iterator it1 = partition->begin(); it1 != partition->end(); ++it1) {
						h.remove_node(*it1);
						current_solution.insert(*it1);
					}
					set<Node> help_solution = run_iter_comp(h);
					current_solution.insert(help_solution.begin(), help_solution.end());
					// found a new best partial solution using the current multi-edge branching nodes
					if (current_solution.size() < partial_solution.size() || partial_solution.size() == 0) {
						partial_solution = current_solution;
					}
				}
			}
			// if there are no multi-edges, just run iterative compression on the connected component
			else {
				debug cout << "There are no multi-edges." << endl;
				partial_solution = run_iter_comp(it);
			}
		}
		partial_solutions.push_back(partial_solution);
		debug cout << "Found partial solution: ";
		debug it.print_nodeset(partial_solutions.back());
		complete_solution.insert(partial_solutions.back().begin(), partial_solutions.back().end());
	}

	debug cout << "----------------------------------------------" << endl;
	debug cout << "Computed a total of " << partial_solutions.size() << " partial solutions with ";
	debug cout << complete_solution.size() << " nodes. Adding " << necessary_nodes.size() << " necessary ";
	debug cout << "nodes we get a solution size of " << (complete_solution.size() + necessary_nodes.size()) << "." << endl;
	debug cout << endl;
	complete_solution.insert(necessary_nodes.begin(), necessary_nodes.end());
	for (const auto &it : complete_solution) {
		cout << node_names.second[it] << endl;
	}
	debug cout << "Sanity check: " << (is_fvs(orig, complete_solution) ? "PASS" : "FAILED") << endl;
	debug cout << "--------------- END OF PROGRAM ---------------" << endl;
	/*
	auto sol2= greedy(orig);
	cout<< "Sol: " << complete_solution.size() << endl;
	cout<< "Sol2: " << sol2.size() << endl;	
	orig.print_nodeset(sol2);
	if(sol2.size()<complete_solution.size())
	{
		while(true){}
	}*/
}

