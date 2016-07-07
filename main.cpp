
#include "fvs_solver.hpp"

using namespace fvs;
using namespace FvsGraph;
#ifndef debug
    #define debug if(true)
#endif
#define MAX_OUTPUT_ARTICULATION 40

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
	pair<unordered_set<Node>, unordered_set<Edge> > master_of_arts = g.get_articulation_elements();
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
	debug cout << "----------------------------------------------" << endl;
	debug cout << "Starting iterative compression" <<endl;
	debug cout << "----------------------------------------------" << endl;
	// "preprocessing"
	g.delete_low_degree_nodes();
	// start! :)
	set<Node> min_fvs = compute_min_fvs(g);
	// sanity check and output of the results
	if (is_fvs(g, min_fvs)) {
		debug cout << "[Partial] It did find a minimal FVS of size " << min_fvs.size() << ". " << endl;
		//print_nodes(min_fvs);
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

void contract_edges(Graph &g) {
    Node u = INVALID_NODE;
    Node v = INVALID_NODE;
    int contracted_edges = 0;
    //debug cout << "Contracting edges: ";
    for(auto &it : g.get_low_degree_nodes()) {
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
            
            if(!g.has_edge(u, v)) {
                //debug cout << "("<<g.get_node_name(u)<<"<->"<<g.get_node_name(it)<<"<->"<<g.get_node_name(v)<<") => ("<<g.get_node_name(u)<<"<->"<<g.get_node_name(v)<<"),";
                ++contracted_edges;
                g.add_edge(u, v);
                g.remove_node(it);
            }
             
        }
    }
    //debug cout <<endl<<endl;
    debug cout << "Contracted " << contracted_edges << " edges." <<endl;
}

void get_connected_graphs(const Graph &g, const list<set<Edge> > &connected_components, list<Graph> &connected_graphs) {
    int i = 0;
    for(const auto &cc : connected_components) {
        Graph* h = new Graph();
        /*cout << "Adding: "<<endl;

        for(const auto &ccc : cc) {
            ++i;
            //cout << "("<<g.get_node_name(ccc.first) <<","<<g.get_node_name(ccc.second)<<"),";	
        }
        cout <<endl;*/
        h->add_edges(cc);
        h->assign_names(g.get_mapping());
        
        i+= cc.size(); 
        // Confused that this number seems unreasonably high?
        // Keep in mind, that we ignore double edges, which are counted multiple times here.
        // Other than that, paths are ignored, because trivial.
        
        if(h->get_n() > 0 && h->get_m() > 0)
        connected_graphs.push_back(*h);
    }
    
    
  	debug cout << "Found/Created "<< connected_components.size() << " connected components with "<<i<<" edges in total." <<endl;
}

int main(int argc, char** argv) {
    // Read graph and store information in variables.
    
    GraphData graph_data = read_graph();
    
    set<Node> necessary_nodes = graph_data.necessary_nodes;
    Mapping node_names = graph_data.mapping;
    Graph g = graph_data.graph;
    
    g.assign_names(node_names);
    
    // Create original copy.
    Graph orig(g);

    // preprocessing
    remove_bridges(g, g.get_articulation_elements().second);
    contract_edges(g);

    // iteratively remove semidisjoint cycles
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
    for(auto &it : connected_graphs) {
	    debug cout << "Trying to find partial solution for graph "<<it.get_name()<<" [n="<<it.get_n()<<"|m="<<it.get_m()<<"]"<<endl;	 
        
        partial_solutions.push_back(run_iter_comp(it));
        debug cout << "Found partial solution: ";
        debug it.print_nodeset(partial_solutions.back());
        complete_solution.insert(partial_solutions.back().begin(), partial_solutions.back().end());
    }
    
    debug cout << "----------------------------------------------"<<endl;
    debug cout << "Computed a total of "<<partial_solutions.size()<<" partial solutions with ";
    debug cout << complete_solution.size() << " nodes. Adding " << necessary_nodes.size() << " necessary ";
    debug cout << "nodes we get a solution size of " << (complete_solution.size()+necessary_nodes.size())<<"."<<endl;
    debug cout << endl;
    complete_solution.insert(necessary_nodes.begin(), necessary_nodes.end());
    for(const auto &it : complete_solution) {
        cout << node_names.second[it] << endl;   
    }
    debug cout << "Sanity check: " << (is_fvs(orig, complete_solution)?"PASS":"FAILED")<<endl;		
    debug cout << "--------------- END OF PROGRAM ---------------" << endl;
}
