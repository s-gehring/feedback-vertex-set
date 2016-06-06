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
#define SKIP_MULTIEDGES
#define TEST_NUMBER 2000000



namespace FvsGraph {

  typedef int Node;
  typedef std::pair<int,int> Edge;
  typedef std::unordered_map<Node, bool> Neighborhood; 
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

          int multiedges;
          std::set<Edge> all_multiedges;
          Sizes sizes;
          
          std::string toString2(const NodeVector &x) {
              std::string res = "";
              bool first = true;
              int s = 0;
              for(NodeVector::const_iterator it = x.begin(); it != x.end(); ++it ) {
                  ++s;
                  if(first) {
                    first = false;
                    res += std::to_string(*it);
                  } else {
                    res += ","+std::to_string(*it);
                  }
                  if(s == 30) return res + std::string(", ...");
              }
              return res;
          }
          std::string toString(const Neighborhood &x) {
              std::string res = "";
              bool first = true;
              for(Neighborhood::const_iterator it = x.begin(); it != x.end(); ++it ) {
                  if(first) {
                    first = false;
                    res = std::to_string(it->first)+std::string(it->second?"M":"S");
                  } else {
                    res += std::string(",")+std::to_string(it->first)+std::string(it->second?"M":"S");
                  }
              }
              return res;
          }
          
          
      public:
          
          Graph() {
            n = m = multiedges = 0;
          }
          
          void clear_node(const Node v) {
            std::set<Node> targets;
            for(const auto& it : adj[v]) {
              targets.insert(it.first);
            }
            for(const auto& it : targets) {
              remove_edge(v, it);
            }
          }
          
          bool has_node(const Node v) { return adj.find(v) != adj.end(); }                                        // O(1)
          bool has_edge(const Node u, const Node v) { return (has_node(u)) && (adj[u].find(v) != adj[u].end()); } // O(1) + O(1)
          bool has_multiedge(const Node u, const Node v) { return has_edge(u, v) && adj[u][v]; }                  // O(1) + O(1)
          bool has_multiedge(const Edge &e) { return has_multiedge(e.first, e.second); }                          // Alias to (const Node, const Node)
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
                      Node w = neighbor->first;
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
              if(multiedges > 0) return true;
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
                      Node w = neighbor->first;
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
              if(m < n) return std::pair<std::list<Node>, bool>(std::list<Node>(), false);
              #ifndef SKIP_MULTIEDGES
              if(multiedges > 0) {
                // Return this multiedge. Nothing fancy.
                std::list<Node> x;
                Edge e = *(all_multiedges.begin());
                x.push_back(e.first);
                x.push_back(e.second);
                return std::pair<std::list<Node>, bool>(x, true);
              }
              #endif
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<std::pair<Node, Node> > S;
           
	bool has_cycle(Graph& g);   for(AdjacencyList::iterator it = adj.begin(); it != adj.end(); ++it) {
                if(done.find(it->first) == done.end()) {
                  S.push(std::make_pair(it->first, INVALID_NODE));
                  done.insert(it->first); 
                  while(!S.empty()) {
                    Node v = S.top().first;
                    Node f = S.top().second;
                    S.pop();
                    for(Neighborhood::const_iterator neighbor = adj[v].begin(); neighbor != adj[v].end(); ++neighbor) {
                      Node w = neighbor->first;
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
                if(it->second) {
                    std::set<Edge>::iterator e = all_multiedges.find(Edge(u, it->first));
                    if(e != all_multiedges.end()) {
                        --multiedges;
                        --m;
                        all_multiedges.erase(e);
                    }
                }
                adj[it->first].erase(u);
                --m;
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
            int added = 0;
            std::unordered_set<Node> s;
            for(std::list<Edge>::const_iterator it = E.begin(); it != E.end(); ++it) {
                add_node(it->first);
                add_node(it->second);
                if(it->first == it->second) continue;
                if(has_multiedge(it->first, it->second)) continue;
                s.insert(it->first);
                s.insert(it->second);
                ++added;
                ++m;
                if(has_edge(it->first, it->second)) {
                    adj[it->first][it->second] = true;
                    adj[it->second][it->first] = true;
                } else {
                    adj[it->first][it->second] = false;
                    adj[it->second][it->first] = false;
                    all_multiedges.insert(*it);
                    ++multiedges;
                }
                
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
                //increment_size(u);
                //increment_size(v);
                adj[u][v] = false;
                adj[v][u] = false;
              } else {
                if(has_multiedge(u,v)) {
                  return false;
                }

                //increment_size(u);
                //increment_size(v);
                adj[u][v] = true;
                adj[v][u] = true;
                all_multiedges.insert(std::pair<Node, Node>(u,v));
                ++multiedges;
              }

              ++m;

              
              sizes.insert(adj[u]);
              sizes.insert(adj[v]);
               
              return true;
          }
          
          void clear() {
            adj.clear();
            sizes.clear();
            multiedges = n = m = 0;
            all_multiedges.clear();
            sizes.clear();
          }
          
          bool remove_edge(const Node u, const Node v) {
            if(!has_node(u) || !has_node(v)) return false;
            if(adj[u].find(v) == adj[u].end()) return false;
            sizes.erase(adj[u]);sizes.erase(adj[v]);
            adj[u].erase(v);
            adj[v].erase(u);
            --m;
            std::set<Edge>::iterator e = all_multiedges.find(Edge(u, v));
            if(e != all_multiedges.end()) {
                --multiedges;
                --m;
                all_multiedges.erase(e);
            }
            sizes.insert(adj[u]);sizes.insert(adj[v]);
            return true;
          }

          
          
          void print_all_edges() {
              for(AdjacencyList::iterator it = adj.begin(); it != adj.end(); ++it) {
                  // it-> first = key;
                  // it-> second= value;
                  std::string s = toString(it->second);
                  printf("[%i:[%s]]\n", it->first, s.c_str());
              }
          }
  };
  
}

#endif
