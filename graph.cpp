#ifndef _GRAPH_CPP
#define _GRAPH_CPP

#include <vector>
#include <time.h>

#include <stack>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <map>
#include <string>
#include <random>
#include <math.h>

#include "graph.hpp"




namespace FvsGraph{

  
          bool Graph::is_deg_three() const {
            for(const auto &it :get_low_degree_nodes()) {
              if(get_single_degree(it) < 3) return false; 
            }
            return true;
          }
  
          int Graph::get_n() const {
            return n;
          }
          int Graph::get_m() const {
            return m; 
          }
          

          std::string Graph::get_name() const {
            char s[80];
            sprintf(s, "%p", (void*) this);
            return std::string(s);
          }
    
  
          Graph::Graph() {
            #ifdef __DEBUG
            d = Debugger::get_instance("logs/graph.log", Debugger::ALL);
            d->log("Created graph with pointer " + get_name() + ".", Debugger::DEBUG);
            #endif
            n = m;
            
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
            #ifdef __DEBUG
            debug("Deleted "+std::to_string(changes)+ " nodes with degree at most one.");
            #endif
            return changes;
          }
            
  
          
          void Graph::clear_node(const Node v) {
            std::set<Node> targets;
            std::string tmp = "Clearing node "+std::to_string(v)+ ": ";
            if(!has_node(v)) {
              #ifdef __DEBUG
              warn(tmp+ "Node doesn't exist."); 
              #endif
            } else {
              for(const auto& it : adj[v]) {
                targets.insert(it);
              }
              for(const auto& it : targets) {
                remove_edge(v, it);
              }
              #ifdef __DEBUG
              note(tmp+"Cleared "+std::to_string(targets.size())+" edges.");
              #endif
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
                  #ifdef __DEBUG
                  warn("Trying to find a semidisjoint cycle, and found a node with degree 1: "+std::to_string(current_node)+" this will still work, but may need some unnecessary time.");
                  #endif
                  not_done = false;
                } else if(get_single_degree(current_node) > 2) {
                  not_done = false;
                } else {
                  if(current_node == semi_disjoint_path_one.back()) {
                    #ifdef __DEBUG
                    warn("Found a full disjoint cycle. I don't really believe this.");
                    #endif
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
                  #ifdef __DEBUG
                  warn("Trying to find a semidisjoint cycle, and found a node with degree 1: "+std::to_string(current_node)+" this will still work, but may need some unnecessary time.");
                  #endif
                  not_done = false;
                } else if(get_single_degree(current_node) > 2) {
                  not_done = false;
                } else {
                  if(current_node == semi_disjoint_path_two.front()) {
                    // This is not possible.
                    #ifdef __DEBUG
                    err("Logic error. Got to source node without returning before.");
                    #endif
                    return std::make_pair(semi_disjoint_path_two, true);
                  }
                }
              }
              
              if(semi_disjoint_path_one.front() == semi_disjoint_path_two.back()) {
                semi_disjoint_path_two.pop_back();
				        semi_disjoint_path_two.pop_front();
                semi_disjoint_path_one.splice(semi_disjoint_path_one.end(), semi_disjoint_path_two);
                #ifdef __DEBUG
                debug("Found a disjoint cycle of size "+std::to_string(semi_disjoint_path_one.size())+".");
                #endif
                return std::make_pair(semi_disjoint_path_one, true);
              }
            }
            return std::make_pair(std::list<Node>(), false);
          }
          
          bool Graph::has_cycle() const {            
              if(m >= n) {
                #ifdef __DEBUG
                note("Has_Cycle: Use heuristic: There is a cycle, because m [="+std::to_string(m)+"] >= n [="+std::to_string(n)+"]");
                #endif
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
                          #ifdef __DEBUG
                          note(tmp+"found a cycle.");
                          #endif
                          return true;
                        }
                      }
                    }
                  }
                }
              }
              #ifdef __DEBUG
              note(tmp+" found no cycle.");
              #endif
              return false;
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
                          bool not_done = true;
                          while(not_done) {
                            if(pre.find(v) == pre.end()) {
                              #ifdef __DEBUG
                              err("Found a cycle, but can't follow back to original node.");
                              #endif
                              return std::pair<std::list<Node>, bool>(std::list<Node>(), true);
                            }
                            v = pre[v];
                            if(v == INVALID_NODE) {
                              err("Found a cycle, but can't follow back to original node. Found invalid node."); 
                              return std::pair<std::list<Node>, bool>(std::list<Node>(), true);
                            }
                            if(v == w) {
                               #ifdef __DEBUG
                               note("Found a cycle, and found its origin. Returning cycle.");
                               #endif
                               x.push_back(w);
                               not_done = false;
                               return std::pair<std::list<Node>, bool>(x, true);
                            }
                            x.push_back(v);
                          }
                          #ifdef __DEBUG
                          note(tmp+"found a cycle.");
                          #endif
                          return std::pair<std::list<Node>, bool>(x, true);
                        }
                      }
                    }
                  }
                } // else { Node has already been found. }
              }
            #ifdef __DEBUG
              note(tmp+"found no cycle.");
            #endif
              return std::make_pair<std::list<Node>, bool>(std::list<Node>(), false);
          }
          
          int Graph::articulate(const Node u, bool vis[], int dsc[], int low[], int par[], std::unordered_set<Node> &a_n, std::unordered_set<Edge> &a_e, int time) const {
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
          
          std::pair<std::unordered_set<Node>, std::unordered_set<Edge> > Graph::get_articulation_elements() const {
            // TODO: n = max|V| over all time
            bool *vis = new bool[n];
            int *dsc = new int[n];
            int *low = new int[n];
            int *par = new int[n];
            
            std::unordered_set<Node> art_nodes;
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
            if(candidates.size() == 0) {
              #ifdef __DEBUG
              warn("lowest_deg_node called with an empty set. Returning invalid node (-1).");
              #endif
            }
            // Iterate over candidates
            Node minCan = *(candidates.begin());
            size_t min = adj.find(minCan)->second.size();
            
            for(const auto& it : candidates) {
                if(adj.find(it)->second.size() < min) {
                    min = adj.find(it)->second.size();
                    minCan = it;
                }
            }
            return minCan;
          }
          
          bool Graph::add_node(const Node u) { // O(1)
            if(has_node(u)) {
              #ifdef __DEBUG
              note("Added node which is already existant: "+std::to_string(u));
              #endif
              return false;
            }
            adj[u] = Neighborhood();
            ++n;
            low_deg_nodes.insert(u);
            #ifdef __DEBUG
            debug("Added note "+std::to_string(u));
            #endif
            return true;
          }
          
          std::pair<Neighborhood, bool> Graph::get_neighbors(const Node u) const {
          /* O(1) */
              if(!has_node(u)) {
                #ifdef __DEBUG
                warn("Requested Neighborhood of a node, which doesn't exist. Node("+std::to_string(u)+")");
                #endif
                return std::pair<Neighborhood, bool>(Neighborhood(), false);
              }
              #ifdef __DEBUG
              debug("Requested Neighborhood of node "+std::to_string(u)+". Node has degree "+std::to_string(get_single_degree(u))+".");
              #endif
              return std::pair<Neighborhood, bool>(adj.find(u)->second, true);
          }
          
          bool Graph::remove_node(const Node u) {
            if(!has_node(u)) {
              #ifdef __DEBUG
              warn("Removing node, which doesn't exist: Node("+std::to_string(u)+").");
              #endif
              return false;
            }
            low_deg_nodes.erase(u);
            // Fuck. Go to each neighbor and inform him about the change.
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
              #ifdef __DEBUG
              err("Get degree of a non-existant node "+std::to_string(u)+".");
              #endif
              return -1;
            }
            return adj.find(u)->second.size();
          }
  
          
          void Graph::add_edges(const std::list<Edge> &E) {
            std::unordered_set<Node> s;
            for(const auto &it : E) {
                add_node(it.first);
                add_node(it.second);
                if(it.first == it.second) continue;
                
                s.insert(it.first);
                s.insert(it.second);
                
                ++m;
                
                adj[it.first].insert(it.second);
                adj[it.second].insert(it.first);
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
                #ifdef __DEBUG
                warn("Trying to add recursive edge, which is not allowed. Recursion on node "+std::to_string(u)+".");
                #endif
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
                #ifdef __DEBUG
                warn("Trying to add existant edge. Edge in question: ("+std::to_string(u)+"|"+std::to_string(v)+").");
                #endif
                return false;
              }
          }
          
          void Graph::clear() { 
            adj.clear();
            n = m = 0;
            low_deg_nodes.clear();
            #ifdef __DEBUG
            note("Clearing Graph.");
            #endif
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
              #ifdef __DEBUG
              warn("Trying to remove edge ("+std::to_string(u)+"|"+std::to_string(v)+"), but at least one of its nodes doesn't exist.");
              #endif
              return false;
            }
            if(adj[u].find(v) == adj[u].end()) {
              #ifdef __DEBUG
              warn("Trying to remove edge ("+std::to_string(u)+"|"+std::to_string(v)+"), but it doesn't exist.");
              #endif
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
