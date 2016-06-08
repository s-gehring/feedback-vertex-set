#ifndef _GRAPH_CPP
#define _GRAPH_CPP

#include <boost/container/flat_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
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





#define INVALID_NODE -1



namespace FvsGraph {

  typedef int Node;
  typedef std::pair<int,int> Edge;
  typedef std::unordered_set<Node> Neighborhood; 
  typedef std::unordered_map<Node, Neighborhood> AdjacencyList;  // A hashtable.
  
  typedef std::unordered_set<Node> NodeVector;
  struct compareNeighborhoods {
      bool operator()(const Neighborhood &u,const Neighborhood &v) {
        return u.size() < v.size();
      };
    };
  typedef std::multiset<Neighborhood, compareNeighborhoods> Sizes;  


  class Graph {
      public:
          int n;
          int m;

      private:
          AdjacencyList adj;          

          Sizes sizes;
          
          
      public:
          
          Graph() {
            n = m;
          }
          
          void clear_node(const Node v) {
            std::set<Node> targets;
            for(const auto& it : adj[v]) {
              targets.insert(it);
            }
            for(const auto& it : targets) {
              remove_edge(v, it);
            }
          }
          
          bool has_node(const Node v) { return adj.find(v) != adj.end(); }                                        // O(1)
          bool has_edge(const Node u, const Node v) { return (has_node(u)) && (adj[u].find(v) != adj[u].end()); } // O(1) + O(1)
          bool has_edge(const Edge &e) { return has_edge(e.first, e.second); }                                    // Alias to (const Node, const Node)
          AdjacencyList get_adjacency_list() {
            return adj;
          }
          
          
          std::pair<std::set<Node>, bool> find_semidisjoint_cycle() {
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<std::pair<Node, Node> > S;
              std::set<Node> sd_cycle;
              std::vector<Edge> foundEdges;
              for(AdjacencyList::iterator it = adj.begin(); it != adj.end(); ++it) {
                if(done.find(it->first) == done.end()) {
                  // connected_components++;
                  S.push(std::make_pair(it->first, INVALID_NODE));
                  done.insert(it->first); 
                  while(!S.empty()) {
                    Node v = S.top().first;
                    Node f = S.top().second;
                    S.pop();
                    for(Neighborhood::const_iterator neighbor = adj[v].begin(); neighbor != adj[v].end(); ++neighbor) {
                      Node w = *neighbor;
                      if(done.find(w) == done.end()) {
                        // Tree edge
                        foundEdges.push_back(std::make_pair(v, w));
                        S.push(std::make_pair(w, v));
                        done.insert(w);
                      } else {
                        if(w != f) {
                        // Target(e) == w
                        // Source(e) == v
                          if(find(foundEdges.begin(), foundEdges.end(), std::make_pair(w, v)) == foundEdges.end()) {
                            foundEdges.push_back(std::make_pair(v,w));
                            // found cycle, check if semi-disjoint
                            std::vector<std::pair<Node, Node> >::iterator it2 = foundEdges.end() -1;
                            int nodes_with_high_deg = 0;
                            while(it2->first != w) {
                              sd_cycle.insert(it2->second);
                              if(get_single_degree(it2->second) > 2) {
                                nodes_with_high_deg = 0;
                              }
                              --it2;
                            }
                            sd_cycle.insert(it2->second);
                            if(nodes_with_high_deg <= 1) {
                              return std::make_pair(sd_cycle, true);
                            } else {
                              sd_cycle.clear();
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              return std::make_pair(std::set<Node>(), false);
          }
          
          bool has_cycle() {
              if(m >= n) return true;
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
                          return true;
                        }
                      }
                    }
                  }
                }
              }
              return false;
          }
          
          
          std::pair<std::list<Node>, bool> get_cycle() {
              //if(m < n) return std::pair<std::list<Node>, bool>(std::list<Node>(), false);
              
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
                          return std::pair<std::list<Node>, bool>(x, true);
                        }
                      }
                    }
                  }
                }
              }
              return std::make_pair<std::list<Node>, bool>(std::list<Node>(), false);
          }
          
          Node source(const Edge &e) { return e.first; }
          Node src(const Edge &e) {return source(e);}
          Node target(const Edge &e) { return e.second; }
          Node trg(const Edge &e) { return target(e); }
          
          Node lowest_deg_node(const std::set<Node> &candidates) {

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
          
          bool add_node(const Node u) { // O(1)
            if(has_node(u)) return false;
            adj[u] = Neighborhood();
            sizes.insert(adj[u]);
            ++n;
            return true;
          }
          
          std::pair<Neighborhood, bool> get_neighbors(Node u) {
          /* O(1) */
              if(!has_node(u)) return std::pair<Neighborhood, bool>(Neighborhood(), false);
              return std::pair<Neighborhood, bool>(adj[u], true);
          }
          
          bool remove_node(const Node u) {
            if(!has_node(u)) return false;
            sizes.erase(adj[u]);
            
            // Fuck. Go to each neighbor and inform him about the change.
            for(Neighborhood::const_iterator it = adj[u].begin(); it != adj[u].end(); ++it) {
                remove_edge(*it, u);
                // remove_edges()
            }

            adj[u].clear();
            adj.erase(u);
            --n;
            return true;
          }
          
                   
          int get_single_degree(const Node u) {
            return adj[u].size();
          }
          
          void add_edges(const std::list<Edge> &E) {
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
            // Now repair sizes.
            Sizes::iterator found;
            for(std::unordered_set<Node>::const_iterator it = s.begin(); it != s.end(); ++it) {
                found = sizes.find(adj[*it]);
                if(found != sizes.end()) {
                    sizes.erase(found);
                }
                sizes.insert(adj[*it]);
            }
          }
          
          bool add_edge(const Node u, const Node v) {
            /* 
              O(1)
            */
              if(u==v) return false;
              add_node(u);// +1
              add_node(v);// +1
              
              sizes.erase(adj[u]);
              sizes.erase(adj[v]);
              
              if(!has_edge(u,v)) {
                adj[u].insert(v);
                adj[v].insert(u);
              } 

              ++m;

              
              sizes.insert(adj[u]);
              sizes.insert(adj[v]);
               
              return true;
          }
          
          void clear() {
            adj.clear();
            sizes.clear();
            n = m = 0;
            sizes.clear();
          }
          
          bool remove_edge(const Node u, const Node v) {
            if(!has_node(u) || !has_node(v)) return false;
            if(adj[u].find(v) == adj[u].end()) return false;
            
            sizes.erase(adj[u]);
            sizes.erase(adj[v]);
            
            adj[u].erase(v);
            adj[v].erase(u);
            
            --m;
            
            sizes.insert(adj[u]);
            sizes.insert(adj[v]);
            return true;
          }

          
          
          
  };
  
}

#endif
