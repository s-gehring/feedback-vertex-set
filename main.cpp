

#include "fvs_solver.h"

using namespace fvs;

int main(int argc, char** argv) {

    fvs::Graph g;
    fvs::read_graph(g, "mini_graph.txt");
    fvs::print_graph(g);
    
    std::pair<set<fvs::Node>, bool> feedback;
    
    set<fvs::Node> v1, v2;
    typedef boost::graph_traits<fvs::Graph>::vertex_iterator iterator;
    std::pair<iterator, iterator> nIt = boost::vertices(g);
    for (iterator it = nIt.first; it != nIt.second; ++it) {
        v1.insert((*it));
    }
    
    v1 = { 2, 7, 4, 6 };
    v2 = { 0, 1, 3, 5 };
    int k = 4;
    
    Graph h(g);
    feedback = fvs::compute_fvs(h, g, v1, v2, k);
    
    cout << "FVS: " << (feedback.second ? "exists" : "doesnt exists") << endl;
    for (const auto& i : feedback.first) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: " << feedback.first.size() << endl;
    cout << "--------------- END OF PROGRAM ---------------" << endl;
}
