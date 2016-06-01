

#include "fvs_solver.h"

using namespace fvs;

int main(int argc, char** argv) {
    const char* filepath = "mini_graph.txt";
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
    cout << (0.5*v1.size()) << "<= MinFVS"<<endl;
    int k = v1.size();
    int n = num_vertices(g);
    int m = num_edges(g);
    

    int min = k/2;
    int max = n;
    do { // Do binary search??
        k = (min + max) / 2;
        Graph h(g);
        set<fvs::Node> x (v1);
        set<fvs::Node> y (v2);
        feedback = fvs::compute_fvs(h, g, x, y, k);
        // ToFix: compute_fvs seems to destroy some of these values.
        cout << "Finished calculation for k = "<<k<<endl;

        

    } while(!feedback.second && k!=min);
    
    
    
    cout << "FVS: " << (feedback.second ? "exists" : "doesnt exists") << endl;
    for (const auto& i : feedback.first) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: " << feedback.first.size() << endl;
    cout << "--------------- END OF PROGRAM ---------------" << endl;
}
