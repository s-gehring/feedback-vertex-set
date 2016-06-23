#include <utility>
#include <limits.h>
#include <fstream>
#include <sstream>
#include <boost/cstdint.hpp>
#include <boost/variant/variant.hpp>
#include "fvs_solver.hpp"

using namespace fvs;

  void fvs::print_nodes(const set<Node>& s) {
    set<Node>::iterator it = s.begin();
    if (s.size() > 0) {
      cout << "{" << *it;
      while (++it != s.end()) {

        cout << ", " << *it;
      }
      cout << "}" << endl;
    }
    else {
      cout << "{}" << endl;
    }
  }
  
  bool fvs::has_cycle(const Graph& g) { return g.has_cycle(); }

  pair<list<Node>, bool> fvs::find_semidisjoint_cycle(const Graph& g) { return g.find_semidisjoint_cycle(); }

  Node fvs::get_lowest_degree_node(const Graph &g, const set<Node>& u) {
    // Create induced subgraph and find the vertex with lowest degree in there.
    // Only return, if degree is at most one.
    Graph h;
    g.induced_subgraph(h, u);
    for (const auto &it : h.get_adjacency_list()) {
      if (h.get_single_degree(it.first) < 2) {
        return it.first;
      }
    }
    return INVALID_NODE;
  }

  bool fvs::creates_circle(Graph& g, const set<Node>& u, const Node& v) {
    /*
    **  INVALID_NODE doesn't make any cycle.
    */
    if(v == INVALID_NODE) {
      cout << "Warning: Called creates_circle() with invalid node." << endl;
      return false;
    }

    /*
    **  Vertices, which are not in the graph don't make cycles.
    */
    pair<Neighborhood, bool> neighbors = g.get_neighbors(v);
    if(!neighbors.second) {
      cout << "Warning: Called creates_circle() with a vertex (v) not being in the graph (g)."<<endl;
      return false;
    }

    /*
    **  Find the set of all neighbors of the given vertex which are in u.
    */
    set<Node> neighbors_in_u;
    for(const auto &it : neighbors.first) {
        if(u.find(it) != u.end()) {
            neighbors_in_u.insert(it);
        }
    }
    /*
    **  TODO: Finish this!
    **  Currently, this only returns true, if 
    **  v has two neighbors in U, which are also connected.
    **  Either use DFS or induced_subgraph() or both.
    */
    for (const auto& i : neighbors_in_u) {
      pair<Neighborhood, bool> eIt = g.get_neighbors(i); 
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

  Node fvs::two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node> &v) {
    /*
    **  Check every node in U.
    */
    for (const auto& candidate : u) {
      int neighbors = 0;
      /*
      **  Count number of neighbors of such a candidate...
      */
      Neighborhood n = g.get_neighbors(candidate).first;
      for(const auto& j : n) {
          /*
          **  ...which are in V.
          */
          if(v.find(j) != v.end()) {
              if(++neighbors > 1) return candidate;
          }
      }
    }
    return INVALID_NODE;
  }



  pair<set<Node>, bool> fvs::forest_bipartition_fvs(const Graph& orig, Graph& g, set<Node>& f, set<Node>& v2, int k) {
  set<Node> fvs;
  pair<set<Node>, bool> retValue;
  
  if (k < 0 || (k == 0 && has_cycle(g))) {
    //cout << "Returning false, k got small and g has still cycles." << endl;
    return make_pair(fvs, false);
  }

  if (!has_cycle(g)) {
    //cout << "g has no cycle." << endl;
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
      
    //  cout << "Subcall retuned: " << (retValue.second ? "true" : "false") << endl;
    //  cout << "Subcall-Set: "; for (const auto& i : retValue.first) { cout << i << ", "; } cout << endl;
      
      
      if (false == retValue.second) {
    //    cout << "Returning false after subcall. " << endl;
        return make_pair(fvs, false);
      }
      else {
        fvs = retValue.first;
        fvs.insert(w);
        
    //    cout << "Returning from 3.1: ";
        //for (const auto& i : fvs) { cout << i << ", "; } cout << endl;
        return make_pair(fvs, true);
      }
    }
    else {
      Graph h(g);
      f.erase(w);
      h.remove_node(w);
      retValue = forest_bipartition_fvs(orig, h, f, v2, k - 1);

      if(retValue.second) {
        fvs = retValue.first;
        
        fvs.insert(w);
        
      //cout << "Returning: ";
      //  for (const auto& i : f) { cout << i << ", "; } cout << endl;
        return make_pair(fvs, true);
      }
      else {
        v2.insert(w);
        //      cout << "Removing " << w << " from the fvs." << endl;
        return forest_bipartition_fvs(orig, g, f, v2, k);
      }
    }
  }
  else {
    w = get_lowest_degree_node(g, f);
    if (INVALID_NODE != w && orig.get_single_degree(w) < 2) {
//      cout << "Out degree of " << w << " is smaller than 2" << endl;
      Graph h(g);
      f.erase(w);
      h.remove_node(w);

      return retValue = forest_bipartition_fvs(orig, h, f, v2, k);
    } 
    else if (w != INVALID_NODE) {
      f.erase(w);
      v2.insert(w);

//      cout << "Removing " << w << " from the fvs." << endl;
      return forest_bipartition_fvs(orig, g, f, v2, k);
    }
  }
  
  return make_pair(fvs, false);
}

  GraphData fvs::read_graph(const char* filepath) {
    GraphData result;
    int current_node_id = 0;

    /*
    **  Open a file at the given filepath to read from.
    */
    ifstream file(filepath, ios::in);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file.");
    }
    /*
    **  Now, line for line, read the file...
    */
    string line;
    while (getline(file, line)) {
      istringstream iss(line);
      /*
      **  Put the two read strings into src and dst. 
      **  These are hopefully never bigger than 50 Bytes.
      **  If there aren't two strings, seperated by a
      **  space, then throw an error.
      */
      char src[50], dst[50]; 
      int x = sscanf(line.c_str(), "%s %s", src, dst);
      if(x != 2) {
        throw std::runtime_error("Can't parse file. Wrong format?");
      }
      /*
      **  Convert to C++ strings.
      */
      string src_str(src), dst_str(dst);
      
      /*
      **  Find a mapping for the src-Node or
      **  create a new one if it doesn't exist.
      */
      Node s, t;
      if(result.mapping.first.find(src_str) == result.mapping.first.end()) {
        result.mapping.first[src_str] = current_node_id;
        result.mapping.second[current_node_id] = src_str;
        s = current_node_id;
        ++current_node_id;
      } else {
        s = result.mapping.first[src_str];
      }
      
      /*
      **  Find a mapping for the dst-Node or
      **  create a new one if it doesn't exist.
      */
      if(result.mapping.first.find(dst_str) == result.mapping.first.end()) {
        result.mapping.first[dst_str] = current_node_id;
        result.mapping.second[current_node_id] = dst_str;
        t = current_node_id;
        ++current_node_id;
      } else {
        t = result.mapping.first[dst_str];
      }
      /*
      **  We got a reflexive node (edge between a node and itself).
      **  We don't add this edge, since the graph structure doesn't
      **  allow it. However, we save tht.first > *(adj[it.first].begin())e node as necessary (to delete).
      **
      **  If s != t then we add the edge. This will never create multiedges
      **  since the data structure doesn't allow it.
      */
      if(s == t) {
        result.necessary_nodes.insert(s);
      } else {
        result.graph.add_edge(s, t);
      }
    }
    file.close();
    
    /*
    **  If there are nodes, which became necessary after a number of
    **  edge insertions, safely remove all of them from the graph
    **  again.
    */
    for(const auto &it : result.necessary_nodes) {
      result.graph.remove_node(it); 
    }
    
    /*
    **  Return the result object, which is basically a container containing
    **  * graph
    **  * necessary_nodes
    **  * mapping -> which is a pair of two maps (Str->Int and Int->Str)
    */
    return result;
  }

  /*
  **  I deleted this function, since it's already implemented in Graph.delete_low_degree_nodes().
  */
  
  /*
void cleanup(Graph& g)
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
*/

  set<Node> fvs::two_approx_fvs(Graph& orig) {
    Graph g(orig);              // our working copy
    set<Node> f;                 // the approx of the fvs
    map<Node, double> weights;   // individual weights for each node
    stack<Node> s;               // here we save all nodes of the fvs for checking their neccessarity in the end
    g.delete_low_degree_nodes();
    /*
    **  Initialize Weights to always be one
    **  such that nodes with higher degree
    **  are preferred later.
    */
    for (const auto &it : g.get_adjacency_list()) {
      weights[it.first] = 1;
    }
    while (g.get_n() > 0) {
      // contains a semidisjoint cycle?
      // Here, we do not the same as the algorithm in the paper: we only put the node with degree >2 into the fvs
      /*
      **  Check for a semidisjoint cycle 
      **  (a cycle with at most one vertex with degree higher than 2)
      **  and compute it, if such a cycle exists.
      **
      **  This algorithm differs slightly from the corresponding paper
      **  since we pick the node with highest degree.
      **  Notice that there always exists an optimal solution
      **  which picks such a "high-degree" vertex.
      */
      pair<list<Node>, bool> sdcycle = find_semidisjoint_cycle(g);
      if (sdcycle.second) {
        /*
        **  If every node has degree two (ie we have a full
        **  disjoint cycle), then we don't care what
        **  node is going to be picked. Every node will do
        **  the job.
        */
        weights[sdcycle.first.front()] = 0;
        /*
        **  The highest degree node will be in the front of
        **  sdcycle.
        */
        sdcycle.first.pop_front();

        /*
        **  Since we deleted the high degree node,
        **  all the other nodes in the semidisjoint
        **  cycle only create a path.
        **  We can safely delete such a path,
        **  since it will never be part of a cycle.
        */
        for (list<Node>::const_iterator it = sdcycle.first.begin(); it != sdcycle.first.end(); ++it) {
          g.remove_node(*it);
        }
      } else { 
        /*
        **  Now the graph is "clean" from semidisjoint
        **  and full disjoint cycles.
        */
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
      set<Node> to_delete;
      for (const auto &it : g.get_adjacency_list()) {
        if (weights[it.first] < numeric_limits<double>::epsilon()) {
          f.insert(it.first);
          s.push(it.first);
          to_delete.insert(it.first);
        }
      }
      // delete nodes
      for (const auto &it : to_delete) {
        g.remove_node(it);
      }
      /*
      **  Clean up afterwards.
      **  (<=> delete all < 2 degree nodes).
      */
      g.delete_low_degree_nodes();
    }
    /*
    **  Try to shrink the FVS we got.
    **  For this we try to remove each node from
    **  the FVS and see if the resulting set is
    **  still an FVS.
    */
    Node u;
    while (!s.empty()) {
      u = s.top();
      f.erase(u);
      s.pop();
      if (!is_fvs(orig,f)) {
        f.insert(u);
      }
    }
    return f;
  }

  bool fvs::is_fvs(const Graph& g, const set<Node>& fvs) {
    /*
    **  Create a copy, delete all nodes from fvs
    **  and check wether there's a cycle.
    */
    Graph h(g);
    for (set<Node>::iterator it = fvs.begin(); it != fvs.end(); ++it) {
      h.remove_node(*it);
    }
    return !has_cycle(h);
  }

  pair<set<Node>, bool> fvs::compression_fvs(const Graph& orig, const set<Node>& S) {
    Graph g(orig);
    size_t k = S.size() - 1;
    boost::uint_fast64_t n = pow(2, k + 1); // for fvs of large size, this is too small -> need other approach
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
    for (size_t j = 0; j < n; j++) {
      // guess intersection using binary coding
      h = j;
      for (size_t l = 0; l < k + 1; l++) {
        current = h % 2;
        if (current) {
          D.insert(T[l]);
        }
        h = (h - current) / 2;
      }
      // compute G[S\D]
      set<Node> s_without_d = set_minus(S, D);
      Graph g(orig);
      for (const auto &it : g.get_adjacency_list()) {
        if (s_without_d.find(it.first) == s_without_d.end()) {
          g.remove_node(it.first);
        }
      }
      if (!has_cycle(g)) {
        // compute G without D
        Graph g(orig);
        for (set<Node>::iterator it = D.begin(); it != D.end(); ++it) {
          g.remove_node(*it);
        }
        // run forest bipartition fvs
        set<Node> v_without_s = set_minus(V, S);
        result = forest_bipartition_fvs(g,g,v_without_s,s_without_d,k-D.size());
        if (result.second) {
          set<Node> output = set_union(result.first, D);
          return make_pair(output, true);
        }
      }
      D.clear();
    }
    return make_pair(S, false);
  }

  set<Node> fvs::set_minus(const set<Node> S, const set<Node> T) {
    /*
    **  For each item in S...
    */
    set<Node> difference;
    for (set<Node>::iterator it = S.begin(); it != S.end(); ++it){
      /*
      **  ..if it is not in T...
      */
      if (T.find(*it) == T.end()) {
        /*
        **  ...add it to the solution.
        */
        difference.insert(*it);
      }
    }
    return difference;
  }

  set<Node> fvs::set_union(const set<Node> S, const set<Node> T) {
    set<Node> union_;
    union_.insert(S.begin(), S.end());
    union_.insert(T.begin(), T.end());
    return union_;
  }

  set<Node> fvs::compute_min_fvs(const Graph& orig) {
    Graph g(orig);
    // get nodes of the graph
    set<Node> v;
    for (const auto &it : g.get_adjacency_list()) {
      v.insert(it.first);
    }
    // compute 2-approximation
    set<Node> fvs_approx = two_approx_fvs(g);
    cout << "2 approximation of size " << fvs_approx.size() << " is: " << endl;
    orig.print_nodeset(fvs_approx);
      
    pair<set<Node>, bool> fvs_dumb_approx = g.get_dumb_approx(fvs_approx.size());
      
    if(fvs_dumb_approx.second && (fvs_dumb_approx.first.size() < fvs_approx.size())) {
        cout <<"Greedy approximation is better, take it instead!"<<endl;
        fvs_approx = fvs_dumb_approx.first;
    } else {
        cout << "(Greedy approximation is worse with size >= "<<fvs_approx.size()<<")"<<endl;   
    }
      
    // use any subset of half size
    int k = 0.5*(fvs_approx.size()+fvs_approx.size()%2);
    set<Node> v_prime;
    set<Node>::iterator it = fvs_approx.begin();
    for (int i = 0; i < k; i++) {
      v_prime.insert(*it);
      ++it;
    }
    // initialize sets/vectors
    set<Node> v_iter = set_union(v_prime, set_minus(v, fvs_approx)); // =v0
    set<Node> help_nodes = set_minus(fvs_approx, v_prime); // used to create iter_nodes as vector
    vector<Node> iter_nodes;
    for (set<Node>::iterator it = help_nodes.begin(); it != help_nodes.end(); ++it) {
      iter_nodes.push_back(*it);
    }
    pair<set<Node>, bool> result = make_pair(v_prime, false);
    set<Node> f_iter = v_prime; // f_0 = v_prime
    // get the iterative compression going
    for (size_t j = 0; j < fvs_approx.size() - k; j++) {
      f_iter.insert(iter_nodes[j]);
      v_iter.insert(iter_nodes[j]);
      // construct iterative graph
      Graph g(orig);
      set<Node> to_delete = set_minus(v, v_iter);
      for (set<Node>::iterator it1 = to_delete.begin(); it1 != to_delete.end(); ++it1) {
        g.remove_node(*it1);
      }
      // run compression
      result = compression_fvs(g, f_iter);
      if(result.second) {
        f_iter = result.first;
      }
    }
    return f_iter;
  }

  set<Node> fvs::brute_force_fvs(const Graph& orig) {
    Graph g(orig);
    // get nodes of the graph
    vector<Node> v;
    for (const auto &it : g.get_adjacency_list()) {
      v.push_back(it.first);
    }
    // compute 2-approximation
    set<Node> fvs_approx = two_approx_fvs(g);
    cout << "2-approximation of the FVS has size " << fvs_approx.size() << endl;
    int upper_bound = fvs_approx.size();
    int lower_bound = 0.5*(fvs_approx.size() + fvs_approx.size() % 2);
    cout << "The size of the optimal solution must be between " << lower_bound << " and " << upper_bound << "." << endl;
    // start it!
    set<Node> solution = fvs_approx;
    set<Node> guessed_fvs;
    boost::uint_fast64_t num;
    boost::uint_fast64_t t;
    bool current;
    boost::uint_fast64_t h;
    while (lower_bound < upper_bound) {
      bool found = false;
      size_t k = 0.5*(lower_bound + upper_bound); // size of fvs to be considered
      cout << "Considered size: " << k << endl;
      // get first decimal number where k bits are set in binary representation
      num = 0;
      for (size_t j = 0; j < k; j++) {
        num += pow(2, j);
      }
      boost::uint_fast64_t N = pow(2, v.size());
      while (num < N && !found) {
        // guess fvs using binary coding
        h = num;
        for (size_t l = 0; l < v.size(); l++) {
          current = h % 2;
          if (current) {
            guessed_fvs.insert(v[l]);
          }
          h = (h - current) / 2;
          if (h == 0) {
            l = v.size();
          }
        }
        // test if it is an fvs
        if (is_fvs(g, guessed_fvs)) {
          solution = guessed_fvs;
          upper_bound = k;
          found = true;
          cout << "Found FVS of size " << k << ":" << endl;
          print_nodes(solution);
        }
        guessed_fvs.clear();
        // compute next number with k bits set
        // see: https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        t = (num | (num - 1)) + 1;
        num = num == 0 ? 0 : t | ((((t & (~t + 1)) / (num & (~num + 1))) >> 1) - 1);
      }
      if(!found){
        cout << "Did not find an FVS of size " << k << endl;
        if (upper_bound - lower_bound == 2) {
          lower_bound = upper_bound;
        }
        else {
          lower_bound = k;
        }
      }
    }
    return solution;
  }
