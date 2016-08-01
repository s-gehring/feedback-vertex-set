#ifndef _GRAPH_CPP
#define _GRAPH_CPP

#include "graph.hpp"

namespace FvsGraph{
            
           std::set<Edge> Graph::minimal_spanning_forest() const {
            std::unordered_set<Node> done;
            std::stack<Node> s;
            std::set<Edge> mst;
            for(const auto &it : adj) {
                if(done.find(it.first) == done.end()) {
                    s.push(it.first);
                    while(!s.empty()) {
                        Node u = s.top();
                        s.pop();
                        for(const auto &v : get_neighbors(u).first) {
                            if(done.find(v) == done.end()) {
                                if(u > v) {
                                    mst.insert(std::make_pair(v, u));
                                } else {
                                    mst.insert(std::make_pair(u, v));
                                }
                                s.push(v);
                                done.insert(v);
                            }
                        }
                    }
                    
                }
                
            }
              return mst;
          }
    
          bool Graph::reaches(const Node &u, const Node &v) const{
              std::unordered_set<Node> done;
              std::stack<Node> s;
              s.push(u);
              while(!s.empty()) {
                Node u = s.top();
                s.pop();
                for(const auto &x : get_neighbors(u).first) {
                  if(done.find(x) == done.end()) {
                    if(x == v) return true;
                      s.push(x);
                      done.insert(x);
                  }
                }
              }
              return false;
          }
    
          bool Graph::compare_node_degrees(const std::pair<Node, Neighborhood> u, const std::pair<Node, Neighborhood> v) {
               return u.second.size() > v.second.size();
          }
    
          void Graph::print_nodeset(const std::set<Node> U) const {
            bool b = true;
            std::cout << "{";
            for(const auto &it : U) {
                if(b) {
                    std::cout << get_node_name(it);
                    b = false;
                } else {
                    std::cout <<","<<get_node_name(it);
                }
            }
            std::cout <<"}"<<std::endl;
          }
    
          std::list<std::set<Edge> > Graph::get_connected_components() const {
            std::unordered_set<Node> done = std::unordered_set<Node>();
            std::stack<std::pair<Node, Node> > S;

            std::list<std::set<Edge> > result;
            for(const auto &it : adj) {
              if(done.find(it.first) == done.end()) {
                std::set<Edge> cc;
                // Start new connected component.


                S.push(std::make_pair(it.first, INVALID_NODE));
                done.insert(it.first); 
                while(!S.empty()) {
                  Node v = S.top().first;
                  Node f = S.top().second;
                  S.pop();
                  for(const auto &neighbor : adj.find(v)->second) {
                    if(done.find(neighbor) == done.end()) {
                      cc.insert(std::make_pair(neighbor, v));
                      S.push(std::make_pair(neighbor, v));
                      done.insert(neighbor);
                    } else {
                      if(neighbor != f) {
                        cc.insert(std::make_pair(neighbor, v));

                      }
                    }
                  }
                }
                if(cc.size() > 0) 
                  result.push_back(cc);
              }
              
            }
            return result;
          }
  
          bool Graph::assign_names(const Mapping &m) {
            this->mapping = m;   
              return true;
          }
          const Mapping Graph::get_mapping() const {
            return mapping;
          }
    
          bool Graph::is_deg_three() const {
            if(get_low_degree_nodes().size() != get_n()) return false;
            for(const auto &it :get_low_degree_nodes()) {
              if(get_single_degree(it) < 2) return false; 
            }
            return true;
          }
          bool Graph::is_deg_most_three_in_set(const std::set<Node> v1) const {
            for(const auto &u : v1) {
                if (get_neighbors(u).first.size()>3)
                {
                  return false;
                }
            }
              return true;
          }
          
  
          size_t Graph::get_n() const {
            return n;
          }
          size_t Graph::get_m() const {
            return m; 
          }
          

          std::string Graph::get_name() const {
            char s[80];
            sprintf(s, "%p", (void*) this);
            return std::string(s);
          }
    
  
          Graph::Graph() {
            n = m = 0;;
            
          }
  
          std::unordered_set<Node> Graph::get_low_degree_nodes() const {
            return low_deg_nodes;
          }
  
          int Graph::delete_low_degree_nodes() {
            bool changes_occur = true;
            int changes = 0;
            do {
              changes_occur = false;
              std::set<Node> to_remove;
              for(const auto& it : adj) {
                if(it.second.size() < 2) {
                  changes_occur = true;
                  to_remove.insert(it.first);
                }
              }
              changes = changes + to_remove.size();
              for(const auto &it : to_remove) {
                remove_node(it); 
                
              }
               
            } while(changes_occur);
            return changes;
          }
          
          void Graph::clear_node(const Node v) {
            std::set<Node> targets;
            std::string tmp = "Clearing node "+std::to_string(v)+ ": ";
            if(!has_node(v)) {
            } else {
              for(const auto& it : adj[v]) {
                targets.insert(it);
              }
              for(const auto& it : targets) {
                remove_edge(v, it);
              }
            }
          }
          
          bool Graph::has_node(const Node v) const { return adj.find(v) != adj.end(); }                                        // O(1)
          bool Graph::has_edge(const Node u, const Node v) const { return (has_node(u)) && (adj.find(u)->second.find(v) != adj.find(u)->second.end()); } // O(1) + O(1)
          bool Graph::has_edge(const Edge &e) const { return has_edge(e.first, e.second); }                                    // Alias to (const Node, const Node) 
          AdjacencyList Graph::get_adjacency_list() const {
            return adj;
          }
          
          // The DFS algorithms
          std::pair<std::list<Node>, bool> Graph::find_semidisjoint_cycle() const {
            std::unordered_set<Node> no;
            for(const auto &v : low_deg_nodes) { 
              if(no.find(v) != no.end()) continue;
              if(get_single_degree(v) != 2) continue;
              std::list<Node> semi_disjoint_path_one;
              semi_disjoint_path_one.push_back(v);
              no.insert(v);
              Node last_node = INVALID_NODE;
              Node current_node = v;
              
              bool not_done = true;
              while(not_done) {
                Neighborhood::const_iterator it = adj.find(current_node)->second.begin();
                Node next_node = *(it);
                if(next_node == last_node) {
                  std::advance(it, 1);
                  next_node = *(it); 
                }
                last_node = current_node;
                current_node = next_node;
                no.insert(current_node);
                semi_disjoint_path_one.push_front(current_node);
                // TODO: Remove unnecessary check.
                if(get_single_degree(current_node) < 2) {
                  not_done = false;
                } else if(get_single_degree(current_node) > 2) {
                  not_done = false;
                } else {
                  if(current_node == semi_disjoint_path_one.back()) {
					          semi_disjoint_path_one.pop_back();
                    return std::make_pair(semi_disjoint_path_one, true);
                    }
                }
              }
              not_done = true;
              std::list<Node> semi_disjoint_path_two;
			        semi_disjoint_path_two.push_back(v);
              last_node = INVALID_NODE;
              current_node = v;
              while(not_done) {
                Neighborhood::const_iterator it = adj.find(current_node)->second.begin();
                std::advance(it, 1);
                Node next_node = *it;
                if(next_node == last_node) {
                  next_node = *(adj.find(current_node)->second.begin()); 
                }
                last_node = current_node;
                current_node = next_node;
                no.insert(current_node);
                semi_disjoint_path_two.push_back(current_node);
                if(get_single_degree(current_node) < 2) {
                  not_done = false;
                } else if(get_single_degree(current_node) > 2) {
                  not_done = false;
                } else {
                  if(current_node == semi_disjoint_path_two.front()) {
                    // This is not possible.
                    return std::make_pair(semi_disjoint_path_two, true);
                  }
                }
              }
              
              if(semi_disjoint_path_one.front() == semi_disjoint_path_two.back()) {
                semi_disjoint_path_two.pop_back();
				        semi_disjoint_path_two.pop_front();
                semi_disjoint_path_one.splice(semi_disjoint_path_one.end(), semi_disjoint_path_two);
                return std::make_pair(semi_disjoint_path_one, true);
              }
            }
            return std::make_pair(std::list<Node>(), false);
          }
          
          bool Graph::has_cycle() const {            
              if(m >= n && n > 0) {
                return true;
              }
              std::string tmp = "Has_Cycle: Don't use heuristic, because m [="+std::to_string(m)+"] < n [="+std::to_string(n)+"]: DFS ";
              
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<std::pair<Node, Node> > S;
              for(const auto &it : adj) {
                if(done.find(it.first) == done.end()) {
                  S.push(std::make_pair(it.first, INVALID_NODE));
                  done.insert(it.first); 
                  while(!S.empty()) {
                    Node v = S.top().first;
                    Node f = S.top().second;
                    S.pop();
                    for(const auto &neighbor : adj.find(v)->second) {
                      if(done.find(neighbor) == done.end()) {
                        S.push(std::make_pair(neighbor, v));
                        done.insert(neighbor);
                      } else {
                        if(neighbor != f) {
                          return true;
                        }
                      }
                    }
                  }
                }
              }
              return false;
          }
  
          bool Graph::induced_subgraph(Graph &s, const std::set<Node>& u) const {
            s.clear();
            for (const auto& node_to_add : u) {
              std::pair<Neighborhood, bool> eIt = get_neighbors(node_to_add);
              if(!eIt.second) {
                  std::cout << "Warning: Called induced_subgraph with a nodeset containing at least one node not in g"<<std::endl;
              } else {
                for(const auto& j : eIt.first) {
                  if(u.find(j)!=u.end()) {
                    s.add_edge(node_to_add, j);
                  }
                }
              }
            }
            return true;
          }
    
          std::string Graph::get_node_name(const Node u)const {
             if(mapping.second.find(u) != mapping.second.end()) {
                return mapping.second.find(u)->second;
             } else {
                return "["+std::to_string(u)+"]";   
             }
          }
            
          void Graph::print() const {
            std::cout << "Printing a graph ["<<get_name()<<"]" << std::endl;
            std::cout << "Number of nodes: " << get_n() << std::endl;
            std::cout << "Number of edges: " << get_m() << std::endl;
            if(get_m() > 1000 || get_n() > 500) {
              std::cout << "Graph too big, skipping complete printing."<<std::endl;
              return;
            }
            for (const auto &it : get_adjacency_list()) {
              std::cout << "Edges outgoing from " << get_node_name(it.first) << ":" << std::endl;
              for (const auto &eit : it.second) {
                std::cout << get_node_name(it.first) << " -> " << get_node_name(eit) << std::endl;
              }
            }
            std::cout << "---------------------------" << std::endl;
            
          }
  
          void Graph::print_tidy()const {
            std::cout << "Number of nodes: " << get_n() << std::endl;
            std::cout << "Number of edges: " << get_m() << std::endl;
            for (const auto &it : get_adjacency_list()) {
              for (const auto &eit : it.second) {
                if(eit > it.first)
                  std::cout << get_node_name(it.first) << " " << get_node_name(eit) << std::endl;
              }
            }
          }
    
          std::pair<std::list<Node>, bool> Graph::get_cycle() const {
              //if(m < n) return std::pair<std::list<Node>, bool>(std::list<Node>(), false);
              std::string tmp = "Get_Cycle: Don't use heuristic, because not yet implemente. DFS ";
              
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<Node> S;
              std::unordered_map<Node, Node> pre;
              for(const auto &it : adj) {
                if(done.find(it.first) == done.end()) {
                  
                  S.push(it.first);
                  pre[it.first] = INVALID_NODE;
                  done.insert(it.first); 
                  
                  while(!S.empty()) {
                    Node v = S.top();
                    S.pop();
                    for(const auto &neighbor : adj.find(v)->second) {
                      Node w = neighbor;
                      if(done.find(w) == done.end()) {
                        S.push(w);
                        pre[w] = v;
                        done.insert(w);
                      } else {
                        if(w != pre[v]) {
                          // (v, w) creates a cycle.
                          // The cycle lies in v, pre[v], pre[pre[v]], ... pre[...pre[v]...]=w
                          
                          std::list<Node> x;
                        
                          std::unordered_set<Node> used;
                          
                          
                          used.insert(w);
                          used.insert(v);
                          x.push_back(w);
                          x.push_front(v);
                          bool v_done = false;
                          bool w_done = false;
                          bool completely_done = false;
                          while(!completely_done) {
                            if(!v_done) {
                              v = pre[v];
                              if(v == INVALID_NODE || used.find(v) != used.end()) {
                                v_done = true; 
                              } else {
                                x.push_front(v);
                                used.insert(v);
                              }
                            }
                            
                            if(!w_done) {
                              w = pre[w];
                              if(w == INVALID_NODE || used.find(w) != used.end()) {
                                w_done = true; 
                              } else {
                                x.push_back(w);
                                used.insert(w);
                              }
                            }
                            completely_done = v_done && w_done;
                          }
                          // if v is used, then cycle is too long at the back
                          if(used.find(v) != used.end()) {
                            while(v != x.back()) x.pop_back(); 
                          } else
                          if(used.find(w) != used.end()) {
                            while(w != x.front()) x.pop_front(); 
                          }
                          return std::make_pair(x, true);
                        }
                      }
                    }
                  }
                } // else { Node has already been found. }
              }
            
              return std::make_pair<std::list<Node>, bool>(std::list<Node>(), false);
          }
          
         int Graph::articulate(const Node u, std::vector<bool>& vis, std::vector<int>& dsc, std::vector<int>& low, std::vector<int>& par, std::set<Node> &a_n, std::unordered_set<Edge> &a_e, int time) const {
            vis[u] = true;
            dsc[u] = time++;
            int min = dsc[u];
            int children = 0;
            
            
            for(const auto &v : get_neighbors(u).first) {
              if(!vis[v]) {
                ++children;
                par[v] = u;
                int maybe_min = articulate(v, vis, dsc, low, par, a_n, a_e, time);
                
                if(maybe_min < min) min = maybe_min;
                
                
                
                
              } else if(v != par[u]) {
                if(dsc[v] < min) min = dsc[v];
              }
            }
            
            if(min == dsc[u] && par[u] != INVALID_NODE) {
              a_e.insert(std::make_pair(par[u], u));
            }
            
            const bool first  = (par[u] == INVALID_NODE && children > 1);
            const bool second = (par[u] != INVALID_NODE && dsc[par[u]] < min);

            if (first || second) {
              a_n.insert(par[u]);
            }
            
            return min;
          }
          
          std::pair<std::set<Node>, std::unordered_set<Edge> > Graph::get_articulation_elements() const {
            // TODO: n = max|V| over all time
            int n = INVALID_NODE;
            for(const auto &it : adj) {
                if(it.first > n) n = it.first+1;
            }
              
            std::vector<bool> vis(n);
            std::vector<int> dsc(n);
            std::vector<int> low(n);
            std::vector<int> par(n);

            
            std::set<Node> art_nodes;
            std::unordered_set<Edge> art_edges;
            
            for(const auto &i : adj) {
              par[i.first] = INVALID_NODE;
              vis[i.first] = false;
            }
            for(const auto &i : adj) {
              if(!vis[i.first]) {
                articulate(i.first, vis, dsc, low,  par, art_nodes, art_edges, 0);
              }
            }
            return std::make_pair(art_nodes, art_edges);
          }
  
          
  
          Node Graph::source(const Edge &e) const { return e.first; }
          Node Graph::src(const Edge &e) const {return source(e);}
          Node Graph::target(const Edge &e) const { return e.second; }
          Node Graph::trg(const Edge &e) const { return target(e); }
          
          Node Graph::lowest_deg_node(const std::set<Node> &candidates) const {
                Node minCan = INVALID_NODE;
		if(candidates.size() == 0) {
			
		}
		else {
			// Iterate over candidates
			minCan = *(candidates.begin());
			size_t min = adj.find(minCan)->second.size();

			for (const auto& it : candidates) {
				if (adj.find(it)->second.size() < min) {
					min = adj.find(it)->second.size();
					minCan = it;
				}
			}
		}
		return minCan;
          }
          
          bool Graph::add_node(const Node u) { // O(1)
            if(has_node(u)) {
              
              return false;
            }
            adj[u] = Neighborhood();
            ++n;
            low_deg_nodes.insert(u);
            
            return true;
          }
          
          std::pair<Neighborhood, bool> Graph::get_neighbors(const Node u) const {
          /* O(1) */
              if(!has_node(u)) {
                
                return std::pair<Neighborhood, bool>(Neighborhood(), false);
              }
              
              return std::pair<Neighborhood, bool>(adj.find(u)->second, true);
          }
          
          bool Graph::remove_node(const Node u) {
            if(!has_node(u)) {
              return false;
            }
            
            low_deg_nodes.erase(u);
            std::set<Edge> to_remove;
            for(Neighborhood::const_iterator it = adj[u].begin(); it != adj[u].end(); ++it) {
                to_remove.insert(std::make_pair(*it, u));
            }
            remove_edges(to_remove);

            adj[u].clear();
            adj.erase(u);
            --n;
            return true; 
          }
          
          int Graph::get_single_degree(const Node u) const {
            if(!has_node(u)) {
              
              return -1;
            }
            return adj.find(u)->second.size();
          }
  
          
          void Graph::add_edges(const std::set<Edge> &E) {
            std::unordered_set<Node> s;
            for(const auto &it : E) {
                add_node(it.first);
                add_node(it.second);
                if(it.first == it.second) continue;
                
                s.insert(it.first);
                s.insert(it.second);
                if(adj[it.first].insert(it.second).second && adj[it.second].insert(it.first).second) {
                  ++m;
                }
                
            }
            for(const auto &it : s) {
              if(get_single_degree(it) < 4) {
                low_deg_nodes.insert(it); 
              } else {
                low_deg_nodes.erase(it);
              }
            }
          }
          
          bool Graph::add_edge(const Node u, const Node v) {
            /* 
              O(1)
            */
              if(u==v) {
                
                return false;
              }
              if(!has_edge(u,v)) {
                add_node(u);// +1
                add_node(v);// +1


                adj[u].insert(v);
                adj[v].insert(u);

                ++m;
                if(get_single_degree(u) > 3) {
                  low_deg_nodes.erase(u);
                }
                if(get_single_degree(v) > 3) {
                  low_deg_nodes.erase(v); 
                }

                return true;
              } else {
                
                return false;
              }
          }
          
          void Graph::clear() { 
            adj.clear();
            n = m = 0;
            low_deg_nodes.clear();
            
          }
            
          void Graph::remove_edges(const std::set<Edge> &E) {
            std::unordered_set<Node> to_update;
            for(const auto &e : E) {
              if(e.first == e.second) continue;
              if(e.first == INVALID_NODE || e.second == INVALID_NODE) continue;
              if(!has_edge(e)) continue;
              Node u = e.first;
              Node v = e.second;
              
              adj[u].erase(v);
              adj[v].erase(u);
              to_update.insert(u);
              to_update.insert(v);
              --m;
            }
            for(const auto &v : to_update) {
              if(get_single_degree(v) < 4) {
                low_deg_nodes.insert(v); 
              } 
            }
            
          }
    
          bool Graph::remove_edge(const Node u, const Node v) {
            if(!has_node(u) || !has_node(v)) {
              
              return false;
            }
            if(adj[u].find(v) == adj[u].end()) {
              
              return false;
            }
            
            
            adj[u].erase(v);
            adj[v].erase(u);
            
            --m;
            if(get_single_degree(u) < 4) {
              low_deg_nodes.insert(u);
            }
            if(get_single_degree(v) < 4) {
              low_deg_nodes.insert(v); 
            }
            return true;
          }

          

}
#endif
