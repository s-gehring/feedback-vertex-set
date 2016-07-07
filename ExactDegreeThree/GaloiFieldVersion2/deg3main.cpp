#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include "degree3.h"
#include "Tests.h"

int main(int argc, char** argv)
{
	//std::srand(std::time(0));
	GraphData graph_data = fvs::read_graph();
	Graph g= graph_data.graph;
//	Graph g;
//	fvs::read_graph(g, "004.graph");
//	fvs::read_graph(g, "mini_graph.txt");
//	fvs::read_graph(g, argv[1]);
//	fvs::print_graph(g);
	//typedef graph_traits<Graph>::vertex_iterator node_iterator;
	//pair<node_iterator, node_iterator> nIt = vertices(g);
	set<Node> s;
	//for (node_iterator it = nIt.first; it != nIt.second; ++it)
	for (const auto &it : g.get_adjacency_list())
	{
		s.insert(it.first);
	}
	//s.erase(0);
	//s.erase(1);
	Tests test;
	test.testAll();
	solveDegree3(g,s,0);
   /*for(int i=0;i<test.size();i+=2)
   {
	   cout<<"Delete note from "<<source(result.second[i],g) <<endl;
	   remove_vertex(source(result.second[i],g),g);
   }
   print_graph(g);*/
	return 0;
}

