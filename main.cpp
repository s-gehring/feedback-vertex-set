

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
    
    
  
    /*
    **  Test whether the graph has been read correctly.
    */ 
    
    print_graph(g);
    
    /*
    **  Output all articulation vertices and bridges.
    **  Note, that I defined that
    **  "And element is an articulation element iff
    **  it is an articulation node or a bridge."
    **
    **  Basically something we can cut the graph at
    **  to simplify computation (hopefully by
    **  an exponential factor).
    */
    pair<list<Node>, list<Edge> > master_of_arts = g.get_articulation_elements();
    
  
    /*
    **  Output all artiulcation elements.
    */
    cout << "Articulation Vertices:"<<endl <<"{";
    bool first_out = true;
    for(const auto &it : master_of_arts.first) {
      if(first_out) {
        first_out = false;
        cout << it;
      } else {
        cout << ", " << it;
      }
    }
    cout << "}" << endl <<endl;
    
    cout << "Bridges:"<<endl <<"{";
    first_out = true;
    for(const auto &it : master_of_arts.second) {
      if(first_out) {
        first_out = false;
        cout << "("<<it.first<< ","<<it.second<<")";
      } else {
        cout << ", ("<<it.first<< ","<<it.second<<")";
      }
    }
    cout <<"}" <<endl <<endl;
    /*
    **  Do the two approximation, store 
    **  the forest decomposition in v1, v2.
    */
    
    set<Node> v1, v2;
    v1 = two_approx_fvs(g);
    for(const auto &it : g.get_adjacency_list()) {
        // Fill v2
        if(v1.find(it.first)==v1.end()) v2.insert(it.first);
    }
  
    /*
    **  Print out the approximate solution. With some stats.
    */
    cout << "Found approx. solution: "<<endl;
    for(const auto& i : v1) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: "<<v1.size()<<endl;
    cout << (0.5*v1.size()) << " <= MinFVS <= "<< v1.size() <<endl;
    
    /*
    **  Binary search the minimal k in a range
    **  from 'size of approx solution * 0.5' to 'size of approx solution * 1.0'.
    **  TODO: This is the wrong way to do this.
    */ 
    int k = v1.size();
    int min = k/2;
    int max = k;
    std::pair<set<Node>, bool> feedback;
    do {        
        k = (min + max) / 2;
        // TODO: Not copy!
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
    
    
    cout << "Algorithm claims, that he " << (feedback.second ? "found" : "didn't find") << " a feedback vertex set." << endl;
    cout << "However, the sanity function tells us that ";
    if(is_fvs(g, feedback.first)) {
        cout << "it did indeed find a FVS." << endl;
    } else {
        cout << "it didn't find a FVS, because there exists a cycle: ";   
        for(const auto &it : feedback.first) {
            h.remove_node(it);
        }
        list<Node> cycle = h.get_cycle().first;
        cout << "{";
        for(const auto &it : cycle) {
            cout << to_string(it)<< "; ";
        }
        cout << "}"<<endl;
    }
        cout <<endl << "FVS in question: "<<endl;
    for (const auto& i : feedback.first) {
        cout << i << ", ";
    }
    cout << endl;
    cout << "Total size: " << feedback.first.size() << endl;
    cout << "--------------- END OF PROGRAM ---------------" << endl;
}
