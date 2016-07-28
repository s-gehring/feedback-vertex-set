#include "fvs_solver.hpp"

using namespace fvs;
using namespace BinCount;

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

  bool fvs::creates_circle(const Graph& g, const Node v, const vector<int> & nodeToComponent) {
    vector<bool> used (nodeToComponent.size(),false);
    for(const auto &it : g.get_neighbors(v).first) {
      if(nodeToComponent[it] != -1) {
	  	  if (used[nodeToComponent[it]])
	  	  { 
	  	  	return true;
        }
	    	else
	  	  {
	  	  	used[nodeToComponent[it]]=true;
	  	  }
      }
      
    }
    return false;
  }

  Node fvs::two_neighbour_node(const Graph& g, const set<Node> &u, const set<Node> &v) {
    Node result = INVALID_NODE;
    int cur_neighbors = 1;
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
              if(++neighbors > cur_neighbors) {
                  result = candidate;
                  cur_neighbors = neighbors;
              } 
          }
      }
    }
    return result;
  }
  
  void fvs::get_connected_graphs(const Graph &g, const list<set<Edge> > &connected_components, list<Graph> &connected_graphs) {
    int i = 0;
    for(const auto &cc : connected_components) {
        Graph h;
        /*cout << "Adding: "<<endl;
        for(const auto &ccc : cc) {
            ++i;
            //cout << "("<<g.get_node_name(ccc.first) <<","<<g.get_node_name(ccc.second)<<"),";	
        }
        cout <<endl;*/
        h.add_edges(cc);
        h.assign_names(g.get_mapping());
        
        i+= cc.size(); 
        // Confused that this number seems unreasonably high?
        // Keep in mind, that we ignore double edges, which are counted multiple times here.
        // Other than that, paths are ignored, because trivial.
        
        if(h.get_n() > 0 && h.get_m() > 0)
        connected_graphs.push_back(h);
    }
    //debug cout << "Found/Created "<< connected_components.size() << " connected components with "<<i<<" edges in total." <<endl;
  }

  pair<set<Node>, bool> fvs::forest_bipartition_fvs(const Graph& orig, Graph& g, set<Node> v1, set<Node> v2, int k) {
    set<Node> fvs;
    pair<set<Node>, bool> retValue;
    /*
    * no budget but g still has cycles -> return false
    */
    if (k < 0 || (k == 0 && has_cycle(g))) {
      return make_pair(fvs, false);
    }
    /*
    * g is a forest -> return empty set
    */
    if (!has_cycle(g)) {
      return make_pair(fvs, true);
    }
    /*
    * Get partitioning of v2 into connected components.
    */
    stack<Node> s;
    int componentNumber=-1;
    int maxIndex=-1;
    for (const auto &it : g.get_adjacency_list()) {
      if(it.first>maxIndex)
       {
         maxIndex=it.first;
       }
    }
    vector<int> nodeToComponent(maxIndex+1,-1);
    for(const auto &it : v2) {
       if (nodeToComponent[it]==-1)
       {
          componentNumber++;
          nodeToComponent[it]=componentNumber;
          s.push(it);
       }
      while(!s.empty()) {
        Node v = s.top();
        s.pop();
        for(const auto &neigh: g.get_neighbors(v).first) {
          /*
          **  Ignore neighbors, which are not in v2
          */
          if(v2.find(neigh) == v2.end()) continue;
          if(nodeToComponent[neigh] == -1) {
            /*
            **  neigh is in no connected component right now.
            */
            nodeToComponent[neigh]=componentNumber;
            s.push(neigh);
          }
        }
        

      }
    }
    /*
    * Pick a vertex w of v1 which has least two neighbors in v2
    * here, we want to pick the vertex with the highest degree!
    */
	  bool cycle =false;
	  Node w;
    for (const auto& candidate : v1)
    {
    	if (creates_circle(g,candidate,nodeToComponent))
    	{
    		w=candidate;
    		cycle=true;
    		break;
    	}
    }
    if (!cycle)
    {
      if (g.is_deg_most_three_in_set(v1) && degree3)
    	{
        Graph h(g);
        h.delete_low_degree_nodes();
        std::list<Graph> connected_graphs;
        get_connected_graphs(h, h.get_connected_components(), connected_graphs);
        //list<set<Node> > partial_solutions;
        set<Node> complete_solution;
        //int totalFvsSize=0;
        for (auto &it : connected_graphs) {
    		//insert seed
          set<Node> v3;
          for(auto v: v1)
          {
           if (it.get_neighbors(v).second)
           {
             v3.insert(v);
           }
          }
    		  auto subFVS= solveDegree3(it,v3,nodeToComponent);
          complete_solution.insert(subFVS.cbegin(), subFVS.cend());
        }
        degree3=false;
        //auto subFVS2=forest_bipartition_fvs(orig, g,v1, v2,k);
        degree3=true;
        if (fvs.size()+complete_solution.size()<= (unsigned) k)
    		{
    			fvs.insert(complete_solution.cbegin(), complete_solution.cend());
          //if (subFVS2.second!=true)
          {
            //throw;
          }
    			return make_pair(fvs, true);
    		}
    		else
    		{
          //if (subFVS2.second!=false)
          {
            //throw;
          }
    			return make_pair(fvs,false);
    		}
    	}
    	else
    	{
	    	w = two_neighbour_node(g, v1, v2);
	    }
    }
    if (w != INVALID_NODE) {
      /*
      * if both neighbours are in the same connected component, i.e. w creates a cycle in g
      */      
      if (cycle) {
        Graph h(g);
        v1.erase(w);
        h.remove_node(w);
        /*
        * select w and reduce the budget by 1
        */
        retValue = forest_bipartition_fvs(orig, h, v1, v2, k - 1);
        /*
        * if returning NO, then return NO
        */
        if (false == retValue.second) {
          return make_pair(fvs, false);
        }
        /*
        * else add w to the fvs
        */
        else {
          fvs = retValue.first;
          fvs.insert(w);
          return make_pair(fvs, true);
        }
      } else {
        /*
        * both neighbours of w are in different connected components
        */
        Graph h(g);
        v1.erase(w);
        h.remove_node(w);
        /*
        * branch on w
        */
        retValue = forest_bipartition_fvs(orig, h, v1, v2, k - 1);
        if(retValue.second) {
          /*
          * add w to the fvs and reduce the budget by 1
          */
          fvs = retValue.first;
          fvs.insert(w);
          return make_pair(fvs, true);
        }
        else {
          /* 
          * do not select w, move w to v2 -> less connected components in g[v2]
          */
          v2.insert(w);
          return forest_bipartition_fvs(orig, g, v1, v2, k);
        }
      }
    } else {
      /*
      * pick any vertex w that has degree <= 1 in g[v1]
      */
      w = get_lowest_degree_node(g, v1);
      /* 
      * if degree of w <= 1 in the original graph
      * delete w from the current graph and do not select it
      */
      if (INVALID_NODE != w && orig.get_single_degree(w) < 2) {
        Graph h(g);
        v1.erase(w);
        h.remove_node(w);
        return retValue = forest_bipartition_fvs(orig, h, v1, v2, k);
      }
      /*
      * else it has exactly one neighbour in v1 and one in v2
      * -> every cycle going through w does contain its neighbours
      * -> do not branch and move w to v2
      */
      else if (w != INVALID_NODE) {
        v1.erase(w);
        v2.insert(w);
        return forest_bipartition_fvs(orig, g, v1, v2, k);
      }
    }
    return make_pair(fvs, false);
  }

  GraphData fvs::read_graph() {
    GraphData result;
    int current_node_id = 0;

    /*
    **  Open a file at the given filepath to read from.
    */
    /*
    if (!cin.is_open()) {
      throw std::runtime_error("Cannot open file.");
    }
    */
    /*
    **  Now, line by line, read the file...
    */
    string line;
    while (!getline(cin, line).eof()) {
      if(line.empty()) {
          debug cout << "Empty line"<<endl;
          continue;
      }
      if(line.at(0) == '#') {
          debug cout << "Got comment"<<line<<endl;
          continue; // Comments
      }
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
        //std::cout << "src: ["<<src<<"], dst: ["<<dst<<"], line: ["<<line<<"]"<<endl;
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
    **  Create a copy, delete all nodes from the fvs
    **  and check wether there's a cycle.
    */
    Graph h(g);
    for (set<Node>::iterator it = fvs.begin(); it != fvs.end(); ++it) {
      h.remove_node(*it);
    }
    return !has_cycle(h);
  }

  pair<set<Node>, bool> fvs::compression_fvs(const Graph& orig, const set<Node>& S) {
    Graph g_help(orig);
    size_t k = S.size() - 1;
    Bin_count counter(k + 1);
    set<Node> D; // the guessed intersection
    // get nodes of the graph
    set<Node> V;
    for (const auto &it : g_help.get_adjacency_list()) {
      V.insert(it.first);
    }
    // convert set S to a vector
    vector<Node> T;
    for (set<Node>::iterator it = S.begin(); it != S.end(); ++it) {
      T.push_back(*it);
    }
    pair<set<Node>, bool> result;
    while(!counter.is_full()) {
      // guess intersection using the binary counter
      for (size_t l = 0; l < k + 1; l++) {
        if (counter.at(l)) {
          D.insert(T[l]);
	      }
      }
      // compute G[S\D]
      set<Node> s_without_d = set_minus(S, D);
      Graph g1;
      orig.induced_subgraph(g1, s_without_d);
      if (!has_cycle(g1)) {
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
	  counter.increase();
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
    //debug cout << "2 approximation of size " << fvs_approx.size() << " is: " << endl;
    //debug orig.print_nodeset(fvs_approx);
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
      Graph h;
      orig.induced_subgraph(h, v_iter);
      h.delete_low_degree_nodes();
      // run compression
      result = compression_fvs(h, f_iter);
      if(result.second) {
        f_iter = result.first;
      }
    }
    return f_iter;
  }
