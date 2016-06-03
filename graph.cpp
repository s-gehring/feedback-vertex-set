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



std::string B2S(bool b) {
  if(b) return "true";
  return "false";
}
namespace Graph {
  

  typedef int Node;
  typedef std::pair<int,int> Edge;
  typedef std::unordered_map<Node, bool> Neighborhood; 
  typedef std::unordered_map<Node, Neighborhood> SortedVectorVector;  // A hashtable.
    
  
  typedef std::unordered_set<Node> NodeVector;
  struct compareNeighborhoods {
      bool operator()(const Neighborhood &u,const Neighborhood &v) {
        return u.size() < v.size();
      };
    };
  typedef std::multiset<Neighborhood, compareNeighborhoods> Sizes;  


  class Graph {
      public:

      private:
          SortedVectorVector adj;          
          int n;
          int m;
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
          /*
          bool add_size(const Node u, const size_t s) { // O(1)
              if(!has_node(u)) return false;
              size_t l = adj[u].size();
              NodeVector nodes = sizes[l];
              nodes.erase(u); 
              if(sizes.find(l+s) == sizes.end()) {                
                sizes[l+s]=NodeVector();
              }
              sizes[l+s].insert(u);
              return true;        
          }
          
          bool increment_size(const Node u) {
              return add_size(u, 1);
          }
          bool decrement_size(const Node u) {
              return add_size(u,-1);
          }
          */
          /* Implemented iteratively
          bool dfs_visit(
              const SortedVectorVector::iterator &u, 
              std::unordered_map<Node, char> &colors,
              std::unordered_map<Node, Node> &pre,
              std::unordered_map<Node, int> &found,
              std::unordered_map<Node, int> &finish,
              int &time) 
          {
              colors[u->first] = 1;
              ++time;
              found[u->first] = time;
              for(Neighborhood::const_iterator it = u->second.begin(); it != void add_edge_list(const std::List<Edge> &u) {

            for(auto std::List<Edge>::iterator it = u.begin(); it!=u.end(); ++it) {
              add_node(it->first);
              add_node(it->second);
              
              if(!has_multiedge(u, v)) {
                if(!has_edge(u, v)) {
                  adj[u][v] = adj[v][u] = false;
                } else {
                  adj[u][v] = adj[v][u] = true;
                  ++multiedges;
                }
                ++m;
              }
            }
            // Added all edges. Now refresh sizes.
            
          }u->second.end(); ++it ) {
                  if(colors[it->first] == 0) {
                      pre[it->first] = u->first;
                      if(dfs_visit(adj.find(u->first), colors, pre, found, finish, time)) return true;
                  }
                  if(colors[it->first] == 1) return true;
              }
              colors[u->first] = 2;
              ++time;
              finish[u->first] = time;
              return false;
          }
          */
          
          
          
      public:
          Graph(void);
          bool has_node(const Node v) { return adj.find(v) != adj.end(); }                                        // O(1)
          bool has_edge(const Node u, const Node v) { return (has_node(u)) && (adj[u].find(v) != adj[u].end()); } // O(1) + O(1)
          bool has_multiedge(const Node u, const Node v) { return has_edge(u, v) && adj[u][v]; }                  // O(1) + O(1)
          bool has_multiedge(const Edge &e) { return has_multiedge(e.first, e.second); }                          // Alias to (const Node, const Node)
          bool has_edge(const Edge &e) { return has_edge(e.first, e.second); }                                    // Alias to (const Node, const Node)
          
          

          bool has_cycle() {
              if(m >= n) return true;
              if(multiedges > 0) return true;
              std::unordered_set<Node> done = std::unordered_set<Node>();
              std::stack<std::pair<Node, Node> > S;
              for(SortedVectorVector::iterator it = adj.begin(); it != adj.end(); ++it) {
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
              for(SortedVectorVector::iterator it = adj.begin(); it != adj.end(); ++it) {
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
                        all_multiedges.erase(e);
                    }
                }
                adj[it->first].erase(u);
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
              for(SortedVectorVector::iterator it = adj.begin(); it != adj.end(); ++it) {
                  // it-> first = key;
                  // it-> second= value;
                  std::string s = toString(it->second);
                  printf("[%i:[%s]]\n", it->first, s.c_str());
              }
          }
  };
  Graph::Graph(void) {
      n=m=multiedges=0;
      //sizes[0] = NodeVector();
  }

}


void addNewNodes(int count, Graph::Graph &G) {
  for(int i = 0; i < count; ++i) {
    G.add_node(i);
  }
}
void addOldNodes(int count, Graph::Graph &G) {
  for(int i = 0; i < count; ++i) {
    G.add_node(1);
  }
}

void addNewEdges(int count, Graph::Graph &G) {
  int x = sqrt(count);
  std::list<Graph::Edge> e;
  for(int i = 0; i < x; ++i) {
    for(int j = 0; j <= x; ++j) {
      e.push_back(Graph::Edge(i,j));
    }

  }
  G.add_edges(e);
}

void addOldEdges(int count, Graph::Graph &G) {
  for(int i = 0; i < count; ++i) {
    G.add_edge(0, 1);
  }
}

std::string nodeListToString(const std::list<Graph::Node> &x) {
    std::string res = "";
    bool first = true;
    for(std::list<Graph::Node>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if(first) {
            first = false;
            res = std::to_string(*it);
        } else {
            res = res + "|"+ std::to_string(*it);
        }
    }
    return res;
}

bool validate_cycle(std::list<Graph::Node> &x, Graph::Graph &g) {
    for(std::list<Graph::Node>::const_iterator it = x.begin(); it != x.end(); ++it) {
        std::list<Graph::Node>::const_iterator next = std::list<Graph::Node>::const_iterator(it);
        ++next;
        if(next == x.end()) {
            if(!g.has_edge(*it, *(x.begin()))) {
                printf("Is not a cycle, because the following edge is missing: %i, %i\n", *it, *(x.begin()));
                return false;
            }
        } else {
            if(!g.has_edge(*it, *next)) {
                printf("Is not a cycle, because the following edge is missing: %i, %i\n", *it, *next);
                return false;
            }
        }
    }
    return true;
}
/*
int main() {

    Graph::Graph G;
    int i,j;
    int n = TEST_NUMBER;
    time_t s1 = time(NULL);
    printf("Starting adding new nodes.\n");
    addNewNodes(n, G);
    printf("Ending adding new nodes.\n");
    time_t s2 = time(NULL);
    printf("Starting adding old nodes.\n");
    addOldNodes(n, G);
    printf("Ending adding old nodes.\n");
    time_t s3 = time(NULL);
    printf("Starting adding new edges.\n");
    
    addNewEdges(3*n, G);
    
    printf("Ending adding new edges.\n");
    time_t s4 = time(NULL);
    printf("Starting adding old edges.\n");
    addOldEdges(n, G);
    printf("Ending adding old edges.\n");
    time_t s5 = time(NULL);
    printf("Starting checking for cycles.\n");
    bool y = G.has_cycle();
    printf("Ending checking for cycles.\n\n");
    time_t s6 = time(NULL);
    printf("Start getting cycle.\n");
    std::pair<std::list<Graph::Node>, bool> x = G.get_cycle();
    printf("End getting cycle.\n");
    time_t s7 = time(NULL);
    printf("Starting validating cycle.\n");
    y = validate_cycle(x.first, G);
    printf("Ending validating cycle.\n");
    time_t s8 = time(NULL);
    

    printf("n = %i\n", n);
    printf("AddNewNodes:\t%.f\n", difftime(s2,s1));
    printf("AddOldNodes:\t%.f\n", difftime(s3,s2));  
    printf("AddNewEdges:\t%.f\n", difftime(s4,s3));
    printf("AddOldEdges:\t%.f\n", difftime(s5,s4));  
    printf("FindCycle:\t%.f\n", difftime(s6,s5));
    printf("GetCycle:\t%.f\n", difftime(s7,s6));
    printf("ValidateCycle:\t%.f\n", difftime(s8,s7));
    if(x.second) {
        printf("Cycle:\t[%s]\n", nodeListToString(x.first).c_str());
        if(!y) {
            printf("Is not a cycle.\n");
        } else {
            printf("Is a cycle.\n");
        }
    }
    

}
*/



}
