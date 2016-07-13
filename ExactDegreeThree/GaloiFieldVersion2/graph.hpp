
#ifndef _H_GRAPH
#define _H_GRAPH


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
  typedef std::pair<std::map<std::string, Node>, std::map<Node, std::string> > Mapping;
  typedef std::unordered_set<Node> Neighborhood; 
  typedef std::unordered_map<Node, Neighborhood> AdjacencyList;  // A hashtable.

  typedef std::unordered_set<Node> NodeVector;
  struct compareNeighborhoods {
      bool operator()(const Neighborhood &u,const Neighborhood &v) const {
        return u.size() < v.size();
      };
    };
class Graph {
      
      
      public:
          /**
          * @brief Compares tuples of nodes and their neighborhood by their neighborhood size.
          */
          static bool compare_node_degrees(const std::pair<Node, Neighborhood> u, const std::pair<Node, Neighborhood> v);
          
          /**
          * @brief Prints a set of nodes with their original name in a pretty way.
          */
          void print_nodeset(const std::set<Node> U) const;
  
          /**
          * @brief Computes a minimal spanning forest for the current graph.
          *
          * Computes a minimal spanning forest via DFS in linear time.
          * For all edges (u, v) in the resulting set it is guaranteed that
          * u <= v. Keep in mind, that lone vertices (with degree 0) will
          * be ignored in this computation.
          *
          * @returns an edge set, corresponding to the MSF.
          */
          std::set<Edge> minimal_spanning_forest() const;
  
          /**
          * @brief Assign a mapping of names to the nodes of this graph.
          *
          * @returns true
          */
          bool assign_names(const Mapping &m);
  
          /**
          * @brief Computes the induced subgraph in respect to a nodeset.
          *
          * Computes a subgraph in respect to a given nodeset u and
          * saves it in s.
          *
          * @returns true
          */
          bool induced_subgraph(Graph &s, const std::set<Node>& u) const;
  
          /**
          * @brief Returns true iff there exists no node of degree != 3.
          */
          bool is_deg_three() const;
  
          /**
          * @brief Returns true iff there exists no node in v1 of degree > 3.
          */
          bool is_deg_most_three_in_set(const std::set<Node> v1) const;
    
    
          size_t get_n() const;
          size_t get_m() const;
          
          /*
          * @brief Prints the graph in human readable way.
          *
          * Prints the adjacency list of the graph in a (more or less)
          * human readable form. This function uses the given name mapping
          * to translate nodes into string before outputting.
          * However, if such a mapping is not found (or empty for a given
          * node), it will output such "unnamed" nodes in square brackets.
          */
          void print() const;
  
          /*
          * @brief Prints the graph in gephi readable way.
          *
          * Print the adjacenfcy lsit of the graph in a way such that it is
          * easily inserted into a gephi project. However for importing into
          * gephi, keep in mind, that a file needs a "source" and a "target"
          * header, even although we're talking about undirected edges.
          */
          void print_tidy() const;
  
          /*
          * @brief Given a node, returns the assigned name given by an assigned mapping.
          * 
          * If there is an assigned mapping of nodes for this graph, 
          * this function will translate a node (which is basically an 
          * integer) into a string.
          *
          * @returns A string, representing the given node.
          */
          std::string get_node_name(const Node u) const;
    
          /*
          * @brief Computes all connected components of a graph.
          *
          * Returns a list of edgesets, each representing a 
          * connected component on its own. Keep in mind that lone
          * vertices (with degree zero) are not included in this result.
          *
          * @returns A list of edge sets, each corresponding to a connected component.
          */
          std::list<std::set<Edge> > get_connected_components() const;
  
          /**
          * @brief Deletes all nodes with degree at most one.
          *
          * @returns The number of nodes deleted.
          */
          int delete_low_degree_nodes();
  
          /**
          * @brief Returns a unique name for the current graph instance.
          *
          * The name is the pointer to the whole datastructure converted to a string.
          * This name can theoretically be recycled after destroying graphs, but
          * this is extremely unlikely.
          *
          * @returns A name as a string.
          */
          std::string get_name() const;
    
          /**
          * @brief Constructor.
          *
          */
          Graph();
  
          /**
          * @brief Clears node from all adjacent edges.
          *
          * Removes all edges of the given node in time O(|Neighbors|*log|Neighbors|).
          * Can be improved, but doesn't seem necessary.
          *
          * @param [in] v The node to clear.
          * @returns Nothing.
          */
          void clear_node(const Node v);
          
          /**
          * @brief Returns whether v exists in the graph.
          *
          * Returns the existance of a node in time O(1).
          *
          * @param [in] v The node in question.
          * @returns True iff v exists in the graph.
          */
          bool has_node(const Node v) const;
          
          /**
          * @brief Returns whether (u, v) or (v, u) exists in the graph.
          *
          * Returns the existance of an edge in time O(1).
          *
          Sizes sizes;
          * @param [in] u The source node of the edge.
          * @param [in] v The target node of the edge.
          * @returns True iff (u,v) exists in the graph.
          */
          bool has_edge(const Node u, const Node v) const;
          /**
          * @brief Returns whether e exists in the graph.
          *
          * Returns the existance of an edge in time O(1).
          *
          * @param [in] e The edge in question.
          * @returns True iff e exists in the graph.
          */
          bool has_edge(const Edge &e) const;
  
          /**
          * @brief Returns the whole adjacency list of the graph.
          *
          * The adjacency list for iterating. The returned object
          * shouldn't be changed. Runs in O(1).
          *
          * @returns An object of type AdjacencyList
          */
          AdjacencyList get_adjacency_list() const;
              
          /**
          * @brief Returns a semidisjoint cycle if there exists one.
          *
          * Returns a cycle, in which every vertex except for one has
          * degree at most two. For FVS it is easy to see, that
          * you'll only need to include the vertex with highest degree
          * in such a set. However, this function will return false
          * if no such cycle exists. Runtime of O(|V|) if no semi-disjoint
          * cycle exists.std::string Graph::get_node_name(Node u)
          * The front entry in the returned list will be the node with
          * higher degree, if such a node exists.
          *
          * @returns A pair, where the second entry is true iff there is a semidisjoint cycle in the first entry.
          * 
          */
          std::pair<std::list<Node>, bool> find_semidisjoint_cycle() const; 
          
          /**
          * @brief Returns whether the graph has a cycle.
          *
          * Uses heuristics and runs in O(1) most of the time.
          * Worst-Case O(|V|).
          *
          * @returns True iff there exists a cycle.
          */       
          bool has_cycle() const;
  
          /**
          * @brief Returns a cycle in the graph if there exists one.
          *
          * Runs in Worst-Case O(|V|).
          * Can be improved to use heuristics, but doesn't seem necessary.
          *
          * @returns A pair, where the second entry is true iff the first entry is a cycle.
          */       
          std::pair<std::list<Node>, bool> get_cycle() const;
          
          /**
          * @brief Returns the source node of a given edge.
          *
          * @param [in] e The edge.
          * @returns A node
          */       
          Node source(const Edge &e) const;
  
          /**
          * @brief Returns the source node of a given edge.
          *
          * @param [in] e The edge.
          * @returns A node
          */       
          Node src(const Edge &e) const;
  
          /**
          * @brief Returns the target node of a given edge.
          *
          * @param [in] e The edge.
          * @returns A node
          */       
          Node target(const Edge &e) const;
  
          /**
          * @brief Returns the target node of a given edge.
          *
          * @param [in] e The edge.
          * @returns A node
          */       
          Node trg(const Edge &e) const;
          
          /**
          * @brief Returns the node with lowest degree out of a set.
          *
          * This is the naive implementation with runtime
          * O(n). Can be improved, but doesn't seem necessary.
          *
          * @param [in] candidates A set of nodes from which the result is to choose.
          * @returns A node with lowest degree or INVALID_NODE if the set is empty.
          */       
          Node lowest_deg_node(const std::set<Node> &candidates) const;
          
          /**
          * @brief Add a node to the graph.
          *
          * Runs in O(1).
          *
          * @param [in] u The node to add.
          * @returns True iff the graph didn't have this node before.
          */       
          bool add_node(const Node u);
          
          
          /**
          * @brief Gets the whole neighborhood of a node.
          *
          * Runs in O(1). Returns false if the node doesn't exist.
          *
          * @param [in] u The node whose neighborhood to return.
          * @returns A pair whose second entry is true iff the first entry is a neighborhood of u.
          */    
          std::pair<Neighborhood, bool> get_neighbors(const Node u) const;        
  
          
          /**
          * @brief Removes a node from the graph.
          *
          * Runs in O(|Neighbors|). Removes the node and all its
          * adjacent edges.
          *
          * @param [in] u The node to remove.
          * @returns True iff u existed in the graph before.
          */    
          bool remove_node(const Node u);        
  
          /**
          * @brief Returns the degree of a given node.
          *
          * Runs in O(1).
          *
          * @param [in] u The node whose degree to return.
          * @returns The degree of the node if it exists. -1 otherwise.
          */    
          int get_single_degree(const Node u) const;
          
          /**
          * @brief Inserts a list of edges.
          *
          * Runs in O(|list|). However, this function
          * is faster than executing add_edge |list| times.
          *
          * @param [in] E The list of edges to insert.
          * @returns Nothing
          */    
          void add_edges(const std::set<Edge> &E);
          
  
          
          /**
          * @brief Adds a single edge.
          *
          * Runs in O(1) but still slow if executed often.
          * Prefer add_edges().
          *
          * @param [in] u The source node of the edge.
          * @param [in] v The target node of the edge.
          * @returns True iff the edge has been inserted successfully.
          */    
          bool add_edge(const Node u, const Node v);        
  
          /**
          * @brief Clears the whole graph.
          *
          * Runs in O(1).
          *
          */    
          void clear();
            
          /**
          * @brief Removes a set of edges
          *
          * Runs in O(|set|). However, this function
          * is faster than executing remove_edge |set| times.
          *
          * @param [in] E The set of edges to remove.
          * @returns Nothing
          */    
          void remove_edges(const std::set<Edge> &E);
    
          /**
          * @brief Removes a single edge.
          *
          * Runs in O(1) but still slow if executed often.
          * Prefer remove_edges().
          *
          * @param [in] u The source node of the edge.
          * @param [in] v The target node of the edge.
          * @returns True iff the edge has been removed successfully.
          */    
          bool remove_edge(const Node u, const Node v);

          /**
          * @brief Returns bridges and articulation points.
          *
          * Runs in DFS time O(|V|+|E|) and returns
          * two lists. first is the list of all
          * articulation vertices, second is the list
          * of all bridges. Those two lists together form
          * the pair of articulation elements, for which holds:
          *
          * An element of a graph G is an articulation element
          * if and only if removal of this element
          * increases the number of connected components in G.
          */
          std::pair<std::unordered_set<Node>, std::unordered_set<Edge> > get_articulation_elements() const;
  
          
          /*
          * @brief Returns all low degree nodes.
          *
          * Returns an unordered set in constant time, containing
          * all nodes with degree at most three.
          */
          std::unordered_set<Node> get_low_degree_nodes() const;
          const Mapping get_mapping() const;
      private:
          AdjacencyList adj;
          std::unordered_set<Node> low_deg_nodes;
          Mapping mapping;
          int n;
          int m;
          int nodes_with_deg_three;
  
  
          int articulate(const Node u, bool vis[], int dsc[], int low[], int par[], std::unordered_set<Node> &a_n, std::unordered_set<Edge> &a_e, int time) const;
          
  
  
};
class GraphData {
  public:
    Graph graph;
    std::set<Node> necessary_nodes;
    Mapping mapping;
};

}



namespace std {
  template <> struct hash<FvsGraph::Edge> {
    size_t operator()(const FvsGraph::Edge &e) const {
      std::hash<FvsGraph::Node> node_hasher;
      return node_hasher(e.first) ^ node_hasher(e.second);
    }
  }; 
}
#endif
