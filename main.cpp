

#include "fvs_solver.h"

using namespace fvs;

int main(int argc, char** argv) {
    const char* filepath = "graphs/mini_graph.txt";
    if(argc>1) {
        filepath = argv[1];
    }
    fvs::Graph g;
    fvs::read_graph(g, filepath);
    fvs::print_graph(g);
    
    std::pair<set<fvs::Node>, bool> feedback;
    
    set<fvs::Node> v1, v2;
    typedef boost::graph_traits<fvs::Graph>::vertex_iterator iterator;
    std::pair<iterator, iterator> nIt = boost::vertices(g);
    for (iterator it = nIt.first; it != nIt.second; ++it) {
        v1.insert((*it));
    }
    /*
    v1 = { 2, 7, 4, 6 };
    v2 = { 0, 1, 3, 5 };
    */
    v1 = two_approx_fvs(g);
    // Slow, to refactor.
    for(iterator it = nIt.first; it != nIt.second; ++it) {
        if(v1.find(*it)==v1.end()) v2.insert(*it);
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

    do {        
        k = (min + max) / 2;
        // Copy like hell!
        Graph h(g);
        set<fvs::Node> x (v1);
        set<fvs::Node> y (v2);
        
        feedback = fvs::compute_fvs(h, g, x, y, k);
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
    feedback = fvs::compute_fvs(h,g,v1,v2,min);
    
    // TODO: Sanity function.
    
    cout << "Feedback Vertex Set: " << (feedback.second ? "Found" : "Not found!") << endl;
    for (const auto& i : feedback.first) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: " << feedback.first.size() << endl;
    cout << "--------------- END OF PROGRAM ---------------" << endl;
}
