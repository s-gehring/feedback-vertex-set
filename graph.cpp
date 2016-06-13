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

  


          std::string Graph::get_name() {
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
            debug("Deleted "+std::to_string(changes)+ " nodes with degree at most one.");
            
            return changes;
          }
            
  
          
          void Graph::clear_node(const Node v) {
            std::set<Node> targets;
            std::string tmp = "Clearing node "+std::to_string(v)+ ": ";
            if(!has_node(v)) {
              warn(tmp+ "Node doesn't exist."); 
            } else {
              for(const auto& it : adj[v]) {
                targets.insert(it);
              }
              for(const auto& it : targets) {
                remove_edge(v, it);
              }
              note(tmp+"Cleared "+std::to_string(targets.size())+" edges.");
            }
          }
          
          bool Graph::has_node(const Node v) { return adj.find(v) != adj.end(); }                                        // O(1)
          bool Graph::has_edge(const Node u, const Node v) { return (has_node(u)) && (adj[u].find(v) != adj[u].end()); } // O(1) + O(1)
          bool Graph::has_edge(const Edge &e) { return has_edge(e.first, e.second); }                                    // Alias to (const Node, const Node)
          AdjacencyList Graph::get_adjacency_list() {
            return adj;
          }
          
          // The DFS algorithms
          std::pair<std::list<Node>, bool> Graph::find_semidisjoint_cycle() {
            std::unordered_set<Node> no;
            for(const auto &v : adj) { // v.first == Node, v.second == Neighborhood
              if(no.find(v.first) != no.end()) continue;
              if(get_single_degree(v.first) != 2) continue;
              std::list<Node> semi_disjoint_path_one;
              semi_disjoint_path_one.push_back(v.first);
              no.insert(v.first);
              Node last_node = INVALID_NODE;
              Node current_node = v.first;
              
              bool not_done = true;
              while(not_done) {
                Neighborhood::const_iterator it = adj[current_node].begin();
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
                  warn("Trying to find a semidisjoint cycle, and found a node with degree 1: "+std::to_string(current_node)+" this will still work, but may need some unnecessary time.");
                  not_done = false;
                } else if(get_single_degree(current_node) > 2) {
                  not_done = false;
                } else {
                  if(current_node == semi_disjoint_path_one.back()) {
                    warn("Found a full disjoint cycle. I don't really believe this.");
					          semi_disjoint_path_one.pop_back();
                    return std::make_pair(semi_disjoint_path_one, true);
                    }
                }
              }
              not_done = true;
              std::list<Node> semi_disjoint_path_two;
			        semi_disjoint_path_two.push_back(v.first);
              last_node = INVALID_NODE;
              current_node = v.first;
              while(not_done) {
                Neighborhood::const_iterator it = adj[current_node].begin();
                std::advance(it, 1);
                Node next_node = *it;
                if(next_node == last_node) {
                  next_node = *(adj[current_node].begin()); 
                }
                last_node = current_node;
                current_node = next_node;
                no.insert(current_node);
                semi_disjoint_path_two.push_back(current_node);
                if(get_single_degree(current_node) < 2) {
                  warn("Trying to find a semidisjoint cycle, and found a node with degree 1: "+std::to_string(current_node)+" this will still work, but may need some unnecessary time.");
                  not_done = false;
                } else if(get_single_degree(current_node) > 2) {
                  not_done = false;
                } else {
                  if(current_node == semi_disjoint_path_two.front()) {
                    // This is not possible.
                    err("Logic error. Got to source node without returning before.");
                    return std::make_pair(semi_disjoint_path_two, true);
                  }
                }
              }
              
              if(semi_disjoint_path_one.front() == semi_disjoint_path_two.back()) {
                semi_disjoint_path_two.pop_back();
				        semi_disjoint_path_two.pop_front();
                semi_disjoint_path_one.splice(semi_disjoint_path_one.end(), semi_disjoint_path_two);
                debug("Found a disjoint cycle of size "+std::to_string(semi_disjoint_path_one.size())+".");
                return std::make_pair(semi_disjoint_path_one, true);
              }
            }
            return std::make_pair(std::list<Node>(), false);
          }
          
          bool Graph::has_cycle() {            
              if(m >= n) {
                note("Has_Cycle: Use heuristic: There is a cycle, because m [="+std::to_string(m)+"] >= n [="+std::to_string(n)+"]");
                return true;
              }
              std::string tmp = "Has_Cycle: Don't use heuristic, because m [="+std::to_string(m)+"] < n [="+std::to_string(n)+"]: DFS ";
              
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<std::pair<Node, Node> > S;
              for(AdjacencyList::iterator it = adj.begin(); it != adj.end(); ++it) {
                if(done.find(it->first) == done.end()) {
                  S.push(std::make_pair(it->first, INVALID_NODE));
                  done.insert(it->first); 
                  while(!S.empty()) {
                    Node v = S.top().first;
                    Node f = S.top().second;
                    S.pop();
                    for(Neighborhood::const_iterator neighbor = adj[v].begin(); neighbor != adj[v].end(); ++neighbor) {
                      Node w = *neighbor;
                      if(done.find(w) == done.end()) {
                        S.push(std::make_pair(w, v));
                        done.insert(w);
                      } else {
                        if(w != f) {
                          note(tmp+"found a cycle.");
                          return true;
                        }
                      }
                    }
                  }
                }
              }
              note(tmp+" found no cycle.");
              return false;
          }
          
          std::pair<std::list<Node>, bool> Graph::get_cycle() {
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
                    for(const auto &neighbor : adj[v]) {
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
                              err("Found a cycle, but can't follow back to original node.");
                              return std::pair<std::list<Node>, bool>(std::list<Node>(), true);
                            }
                            v = pre[v];
                            if(v == INVALID_NODE) {
                              err("Found a cycle, but can't follow back to original node. Found invalid node."); 
                              return std::pair<std::list<Node>, bool>(std::list<Node>(), true);
                            }
                            if(v == w) {
                               note("Found a cycle, and found its origin. Returning cycle.");
                               x.push_back(w);
                               not_done = false;
                               return std::pair<std::list<Node>, bool>(x, true);
                            }
                            x.push_back(v);
                          }
                          note(tmp+"found a cycle.");
                          return std::pair<std::list<Node>, bool>(x, true);
                        }
                      }
                    }
                  }
                } // else { Node has already been found. }
              }
            
              note(tmp+"found no cycle.");
              return std::make_pair<std::list<Node>, bool>(std::list<Node>(), false);
          }
          
          std::pair<Edge, bool> Graph::get_bridge() {
              return std::make_pair(std::make_pair(INVALID_NODE, INVALID_NODE), false);/*
              std::string tmp = "Get Bridge: Don't use heuristic, because not yet implemente. DFS ";
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<std::pair<Node, Node> > S;
              for(AdjacencyList::iterator it = adj.begin(); it != adj.end(); ++it) {
                if(done.find(it->first) == done.end()) {
                  S.push(std::make_pair(it->first, INVALID_NODE));
                  done.insert(it->first); 
                  while(!S.empty()) {
                    Node v = S.top().first;
                    Node f = S.top().second;
                    S.pop();
                    for(Neighborhood::const_iterator neighbor = adj[v].begin(); neighbor != adj[v].end(); ++neighbor) {
                      Node w = *neighbor;
                      if(done.find(w) == done.end()) {
                        S.push(std::make_pair(w, v));
                        done.insert(w);
                      } else {
                        if(w != f) {
                          S.push(std::make_pair(v,f));
                          std::list<Node> x;
                          while(!S.empty()) {
                            x.push_back(S.top().first);
                            S.pop();
                          }
                          note(tmp+"found a cycle.");
                          return std::pair<std::list<Node>, bool>(x, true);
                        }
                      }
                    }
                  }
                }
              }
            
              note(tmp+"found no cycle.");
              return std::make_pair<std::list<Node>, bool>(std::list<Node>(), false);
              */
          }
  
          // /The DFS algorithms
          Node Graph::source(const Edge &e) { return e.first; }
          Node Graph::src(const Edge &e) {return source(e);}
          Node Graph::target(const Edge &e) { return e.second; }
          Node Graph::trg(const Edge &e) { return target(e); }
          
          Node Graph::lowest_deg_node(const std::set<Node> &candidates) {
            if(candidates.size() == 0) {
              warn("lowest_deg_node called with an empty set. Returning invalid node (-1).");
            }
            // Iterate over candidates
            size_t max = 0;
            Node maxCan = INVALID_NODE;
            for(const auto& it : candidates) {
                if(adj[it].size() > max) {
                    max = adj[it].size();
                    maxCan = it;
                }
            }
            return maxCan;
          }
          
          bool Graph::add_node(const Node u) { // O(1)
            if(has_node(u)) {
              note("Added node which is already existant: "+std::to_string(u));
              return false;
            }
            adj[u] = Neighborhood();
            ++n;
            debug("Added note "+std::to_string(u));
            return true;
          }
          
          std::pair<Neighborhood, bool> Graph::get_neighbors(const Node u) {
          /* O(1) */
              if(!has_node(u)) {
                warn("Requested Neighborhood of a node, which doesn't exist. Node("+std::to_string(u)+")");
                return std::pair<Neighborhood, bool>(Neighborhood(), false);
              }
              debug("Requested Neighborhood of node "+std::to_string(u)+". Node has degree "+std::to_string(get_single_degree(u))+".");
              return std::pair<Neighborhood, bool>(adj[u], true);
          }
          
          bool Graph::remove_node(const Node u) {
            if(!has_node(u)) {
              warn("Removing node, which doesn't exist: Node("+std::to_string(u)+").");
              return false;
            }
            
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
          
          int Graph::get_single_degree(const Node u) {
            if(!has_node(u)) {
              err("Get degree of a non-existant node "+std::to_string(u)+".");
              return -1;
            }
            return adj[u].size();
          }
          
          void Graph::add_edges(const std::list<Edge> &E) {
            std::unordered_set<Node> s;
            for(std::list<Edge>::const_iterator it = E.begin(); it != E.end(); ++it) {
                add_node(it->first);
                add_node(it->second);
                if(it->first == it->second) continue;
                
                s.insert(it->first);
                s.insert(it->second);
                
                ++m;
                
                adj[it->first].insert(it->second);
                adj[it->second].insert(it->first);
                
                
            }
            
          }
          
          bool Graph::add_edge(const Node u, const Node v) {
            /* 
              O(1)
            */
              if(u==v) {
                warn("Trying to add recursive edge, which is not allowed. Recursion on node "+std::to_string(u)+".");
                return false;
              }
              if(!has_edge(u,v)) {
                add_node(u);// +1
                add_node(v);// +1


                adj[u].insert(v);
                adj[v].insert(u);

                ++m;


                return true;
              } else {
                warn("Trying to add existant edge. Edge in question: ("+std::to_string(u)+"|"+std::to_string(v)+").");
                return false;
              }
          }
          
          void Graph::clear() {
            adj.clear();
            n = m = 0;
            note("Clearing Graph.");
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
            
          }
    
          bool Graph::remove_edge(const Node u, const Node v) {
            if(!has_node(u) || !has_node(v)) {
              warn("Trying to remove edge ("+std::to_string(u)+"|"+std::to_string(v)+"), but at least one of its nodes doesn't exist.");
              return false;
            }
            if(adj[u].find(v) == adj[u].end()) {
              warn("Trying to remove edge ("+std::to_string(u)+"|"+std::to_string(v)+"), but it doesn't exist.");
              return false;
            }
            
            
            adj[u].erase(v);
            adj[v].erase(u);
            
            --m;
            
            return true;
          }

          

}
#endif
