
#include "fvs_solver.hpp"

using namespace fvs;
using namespace FvsGraph;
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
void run_iter_comp(const char* filepath) {
	Graph g;
	read_graph(g, filepath);
	// compute min fvs by using iterative compression
	cout << "----------------------------------------------" << endl;
	cout << "Starting iterative compression for file " << filepath << endl;
	cout << "----------------------------------------------" << endl;
	// "preprocessing"
	cleanup(g);
	// start! :)
	set<Node> min_fvs = compute_min_fvs(g);
	// sanity check and output of the results
	if (is_fvs(g, min_fvs)) {
		cout << "It did find a minimal FVS of size " << min_fvs.size() << ". The set consists of: " << endl;
		set<Node>::iterator it = min_fvs.begin();
		cout << *it;
		while(++it != min_fvs.end()) {
			
			cout << ", " << *it;
		}
		cout << endl;
	}
	else {
		cout << "Error: The set we have found is not an FVS!" << endl;
	}
}

/*
* Runs brute force to find a minimum feedback vertex set in a given graph.
*
* @param [in] g The filepath of the file where the graph is in.
*/
void run_brute_force(const char* filepath) {
	Graph g;
	read_graph(g, filepath);
	// compute min fvs by using brute force
	cout << "----------------------------------------------" << endl;
	cout << "Starting brute force for file " << filepath << endl;
	cout << "----------------------------------------------" << endl;
	// "preprocessing"
	cleanup(g);
	// start! :)
	set<Node> min_fvs = brute_force_fvs(g);
	// sanity check and output of the results
	if (is_fvs(g, min_fvs)) {
		cout << "It did find a minimal FVS of size " << min_fvs.size() << ". The set consists of: " << endl;
		set<Node>::iterator it = min_fvs.begin();
		cout << *it;
		while(++it != min_fvs.end()) {
			cout << ", " << *it;
		}
		cout << endl;
	}
	else {
		cout << "Error: The set we have found is not an FVS!" << endl;
	}
}

int main(int argc, char** argv) {
    const char* filepath = "graphs/pace/095.graph";
    if(argc>1) {
        filepath = argv[1];
    }
    Graph g;
    read_graph(g, filepath);
	  //print_graph(g);
	  //output_arti_elems_bridges(g);
	  //run_brute_force(filepath);
		run_iter_comp(filepath);
    cout << "--------------- END OF PROGRAM ---------------" << endl;
}
