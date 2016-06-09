

#include "fvs_solver.h"

using namespace fvs;
using namespace FvsGraph;

int main(int argc, char** argv) {
    const char* filepath = "graphs/mini_graph.txt";
    if(argc>1) {
        filepath = argv[1];
    }
    Graph g;
    read_graph(g, filepath);
    g.clear();
    
    g.add_edge(1,2);
    g.add_edge(2,3);
    g.add_edge(1,3);
    g.add_edge(4,5);
    g.add_edge(5,6);
    g.add_edge(4,6);
    g.add_edge(1,5);
    
    
    
    print_graph(g);
    
    
    std::pair<set<Node>, bool> feedback;
    
    set<Node> v1, v2;

    v1 = two_approx_fvs(g);

    for(const auto &it : g.get_adjacency_list()) {
        if(v1.find(it.first)==v1.end()) v2.insert(it.first);
    }
    cout << "Found approx. solution: "<<endl;
    for(const auto& i : v1) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: "<<v1.size()<<endl;
    cout << (0.5*v1.size()) << " <= MinFVS <= "<< v1.size() <<endl;
    int k = v1.size();
  //  int n = num_vertices(g);
//    int m = num_edges(g);
    

    int min = k/2;
    int max = k;

   // TODO: Do this smarter. This implementation sucks and
    // is not what we want. We want to use the approx solution.
    do {        
        k = (min + max) / 2;
        // Copy like hell!
        Graph h(g);
        set<Node> x (v1);
        set<Node> y (v2);
        
        feedback = fvs::forest_bipartition_fvs(h, g, x, y, k);
        cout << "Finished calculation for k = "<<k<<" [min:"<<min<<"|max:"<<max<<"]"<<endl;
        
        if(feedback.second) {
            max = k;
        } else {
            if(min == k) {
                // Rounding errors;
                k = ++min;
            } else {
                min = k;
            }
            
        }

    } while(max != min);
    // This is a bit weird. Binary search cancels if max == min. This is the size of the min FVS. However, since the loop just breaks,
    // we're not guaranteed a result in currentSolution.
    Graph h(g);
    cout << "Found size of min FVS: "<<min<<", continue to compute min FVS."<<endl;
    feedback = forest_bipartition_fvs(h,g,v1,v2,min);
    
    // Sanity function
    if(is_fvs(g, feedback.first)) {
        cout << "Is a FVS." << endl;
    } else {
        cout << "Is not a FVS." << endl;   
    }
    
    cout << "Feedback Vertex Set: " << (feedback.second ? "Found" : "Not found!") << endl;
    for (const auto& i : feedback.first) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: " << feedback.first.size() << endl;
    cout << "--------------- END OF PROGRAM ---------------" << endl;
}
