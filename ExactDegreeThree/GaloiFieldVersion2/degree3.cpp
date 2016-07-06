#include "degree3.h"

using namespace std;
//using namespace fvs;

int getComponentNumber(const Node & n)
{
	return 0;
}

/*int getID(const set<Node> & u,const  Node & n, vector<int> & nodeID, vector<int> & componentID, int & lastRowUsed)
{
	if (u.find(n) != u.end())
	{
		return nodeID[n];
	}
	else
	{
		if (componentID[getComponentNumber(n)] == -1)
		{
			componentID[getComponentNumber(n)] = lastRowUsed;
			lastRowUsed++;
		}
		return componentID[getComponentNumber(n)];
	}
}*/

/*vector<uint64_t> generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, int row_number, const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow)
{
	vector<uint64_t> result(row_number,0);*/
//	result.zeros();
	/*int mi=min(source(e,g),target(e,g));
	int ma=max(source(e,g),target(e,g));*/
/*	int mi = min(nodeToRow[e.first], nodeToRow[e.second]);
	int ma = max(nodeToRow[e.first], nodeToRow[e.second]);*/
	/*int miN;
	int maN;
	if (nodeToRow[e.first] < nodeToRow[e.second])
	{
		mi = nodeToRow[e.first];
		ma = nodeToRow[e.second];
		miN = e.first;
		maN = e.second;
	}
	else
	{
		ma = nodeToRow[e.first];
		mi= nodeToRow[e.second];
		maN = e.first;
		miN = e.second;
	}*/
	/*edgesUsed[mi][ma]++;
	result[lastVertexIndex[mi][ma]] = 1;

	//cout<<mi<<" "<<ma<<" "<<lastVertexIndex[mi][ma]<<endl;
	if (edgesUsed[mi][ma] == pairNumber[mi] + pairNumber[ma])
	{
		//result(ma,0)=-1;
		result[ma] = 1;
		//cout<<"case 1 "<<ma<<endl;
	}
	else
	{
		//result(lastUsedRow,0)=-1;
		assert(lastUsedRow < row_number);
		result[lastUsedRow] = 1;
		lastVertexIndex[mi][ma] = lastUsedRow;
		lastUsedRow++;
	}
	return result;
}*/


void generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, mat & matrix, int columnNumber ,const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow)
{
	int mi = min(nodeToRow[e.first], nodeToRow[e.second]);
	int ma = max(nodeToRow[e.first], nodeToRow[e.second]);
	edgesUsed[mi][ma]++;
	if (lastVertexIndex[mi][ma] != -1)
	{
		matrix(lastVertexIndex[mi][ma], columnNumber) = 1;
	}
	else
	{
		matrix(mi, columnNumber) = 1;
	}
	if (edgesUsed[mi][ma] == pairNumber[mi] + pairNumber[ma])
	{
		matrix(ma,columnNumber) = 1;
	}
	else
	{
		matrix(lastUsedRow,columnNumber) = 1;
		lastVertexIndex[mi][ma] = lastUsedRow;
		lastUsedRow++;
	}
//	matrix.print("");
}

pair<mat,vector<Edge>> graphToMatrix(const Graph& g, const set<Node>& u)
{
  //typedef graph_traits<Graph>::vertices_size_type id;

  //typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
  set<pair<Edge,Edge>> edgePairs;
  //typedef graph_traits<Graph>::vertex_iterator node_iterator;

  /*map<Node,int> nodeIndex;
  int index=0;*/
 
  /*pair<node_iterator, node_iterator> nIt = vertices(g);
  for (node_iterator it = nIt.first; it != nIt.second; ++it)
  {
    nodeIndex[*it]=index;
    index++;
  }*/

  vector<int> nodeToRow(g.get_n());
  vector<int> componentToRow(g.get_n(), -1);
  int lastUsedRow = 0;
  for (const auto &it : g.get_adjacency_list()) {
	  if (u.find(it.first) != u.end())
	  {
		  nodeToRow[it.first] = lastUsedRow;
		  lastUsedRow++;
	  }
	  else
	  {
		  if (componentToRow[getComponentNumber(it.first)] == -1)
		  {
			  componentToRow[getComponentNumber(it.first)] = lastUsedRow;
			  lastUsedRow++;
		  }
		  nodeToRow[it.first] = componentToRow[getComponentNumber(it.first)];
	  }
  }
  cout <<"Number of Nodes: "<<lastUsedRow << endl;
  vector<int> pairNumber(g.get_n(),0);
  size_t edgeNumber = 0;
  for (const auto& firstNode : u) {
      vector<Edge> neighbours;
      //pair<edge_iterator, edge_iterator> eIt = out_edges(n, g);
      //for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
	  Neighborhood nextToFirstNode = g.get_neighbors(firstNode).first;
	  for (const auto& secondNode : nextToFirstNode) {
		  if (u.find(secondNode)==u.end() || firstNode < secondNode)
		  {
			  edgeNumber++;
		  }
		  Edge e;
		  e.first = firstNode;
		  e.second = secondNode;
		  neighbours.push_back(e);
        //neighbours.push_back(*it);
        //edgePairs.insert(make_pair(nodeIndex[i],nodeIndex[target(*it, g)]));
        //neighbours.push_back(nodeIndex[target(*it, g)]);
      }
      pairNumber[nodeToRow[firstNode]]+=neighbours.size()-1;
      for (int i=0;i<neighbours.size();i++)
      {
        for(int j=0;j<i;j++)
        {
            edgePairs.insert(make_pair(neighbours[i],neighbours[j]));
        }
      }
  }
  cout<<"Number of edges: "<<edgeNumber <<endl << "Number of pairs: " << edgePairs.size() << endl;
  vector<Edge> assignment;
  //int row_number = g.get_n() + 2 * edgePairs.size() - g.get_m();
  vector<vector<int>> edgesUsed;
  vector<vector<int>> lastVertexIndex;
  for (int i=0;i<g.get_n();i++)
  {
    vector<int> v(g.get_n(),0);
//    vector<int> w(g.get_n(),nodeToRow[i]);
	vector<int> w(g.get_n(), -1);
    edgesUsed.push_back(v);
    lastVertexIndex.push_back(w);
  }
  //int lastUsedRow=g.get_n();
  /*for (int i=0;i<num_vertices(g);i++)
  {
    lastVertexIndex[i]=i;
  }*/
  int row_number = lastUsedRow + 2 * edgePairs.size() - edgeNumber;
  cout << "Number of rows: " << row_number << endl<< "Number of columns: "<<2*edgePairs.size()<<endl <<"Number of entries: "<< row_number* 2 * edgePairs.size() <<endl;
  mat matrix(row_number,edgePairs.size()*2);
  int columnNumber=0;
  for(auto & p: edgePairs)
  {
	generateIncidenceVector(g, p.first, edgesUsed, lastUsedRow,matrix ,columnNumber, pairNumber, lastVertexIndex, nodeToRow);
	columnNumber++;
	generateIncidenceVector(g, p.second, edgesUsed, lastUsedRow,matrix, columnNumber, pairNumber, lastVertexIndex, nodeToRow);
	columnNumber++;
    /*auto first=generateIncidenceVector(g,p.first,edgesUsed,lastUsedRow,row_number,pairNumber,lastVertexIndex, nodeToRow);
    auto second=generateIncidenceVector(g,p.second,edgesUsed,lastUsedRow,row_number,pairNumber,lastVertexIndex, nodeToRow);
    join_rows_fast(matrix,first);
    join_rows_fast(matrix,second);*/
    assignment.push_back(p.first);
    assignment.push_back(p.second);
  }
  return make_pair(matrix,assignment);
}

/*vector<int> findRest(const vector <int>& index, int dim)
{
    vector<int> result;
    for(int i=0;i<dim;i++)
    {
      if (find(index.begin(), index.end(), i)==index.end())
      {
        result.push_back(i);
      }
    }
    return result;
}*/

/*mat extractColumns(const mat & input, const vector<int> & index)
{
    mat result;
    for(int i: index)
    {
      result=join_rows(result,input.col(i));
    }
    return result;
}*/

/*vector<int> findOnes(const mat & input, int n)
{
  vector<int> result;
  for (int i=0;i<n;i++)
  {
      for(int j=0;j<input.n_cols;j++)
      {
        if (input.at(i,j)==1)
        {
          result.push_back(j);
          break;
        }
      }
  }
  return result;
}*/

/*mat rearrangeMatrix(const mat & input, const vector<int> & arrangement)
{
  mat result(input.n_rows,input.n_cols);
  for(int i=0;i<arrangement.size();i++)
  {
    result.col(arrangement[i])=input.col(i);
  }
  return result;
}*/


mat transformFullRowRank(mat input)
{
  /*int i=0;
  int inputRank=matRank(input);
  while(inputRank!=input.n_rows)
  {
    mat smallerMat=input;
    smallerMat.shed_row(i);
	int newRank = matRank(smallerMat);
    if (newRank==inputRank)
    {
      input=smallerMat;
    }
    else
    {
      i++;
    }        
  }*/
  input.shed_row(0);
  return input;
}

mat colinearToLinear(const mat & input)
{
  //mat fullRank= transformFullRowRank(input);
  //fullRank.print("transformed row rank");
  mat fullRank = input;
  fullRank.shed_row(0);
  if (fullRank.is_square())
  {
    //return eye<mat>(fullRank.n_rows,fullRank.n_cols);
	  return eye<mat>(0, 0);
  }
  /*mat L, U, P;
  lu(L,U,P,fullRank.t());
  L.print("L");
  U.print("U");
  P.print("P");
  mat res=U.t()*L.t()*P;
  res.print("check");
  vector<int> invSubMatrix= findOnes(P,L.n_cols);*/
  /*vector<int> invSubMatrix;
  mat inverse;
  auto inverseSub= fullRank.inverseSubmatrix();
  std::tie(invSubMatrix, inverse) = inverseSub;*/
  /*vector<int> invSubMatrix =std::get<2>(fullRank.upper_triangle_transform());
  invSubMatrix.resize(fullRank.getHeight());*/
  auto standardForm = fullRank.toStandarForm();
  vector<int> arrangement = standardForm.second;
  vector<int> restIndex;// = findRest(invSubMatrix, fullRank.n_cols);
  for (int i = fullRank.getHeight();i<fullRank.getWidth();i++)
  {
	  restIndex.push_back(i);
  }
  //mat submatrix=fullRank.extractColumns(invSubMatrix);
  //submatrix.print("sub");
  mat newRight = standardForm.first.extractColumns(restIndex);
  //newRight.print("nR");
  //mat inverse=submatrix.i();
//  inverse.print("inverse");
//  mat rest= fullRank.extractColumns(restIndex);
//  rest.print("rest");
//  mat newRight=inverse*rest;
//  newRight.print("newRight");
  //vector<int> arrangement=invSubMatrix;
  //arrangement.insert(arrangement.end(), restIndex.begin(), restIndex.end() );
  mat finalMatrix= join_rows(newRight.t(),eye<mat>(fullRank.n_cols-fullRank.n_rows, fullRank.n_cols-fullRank.n_rows));
  //finalMatrix.print("finalMatrix");
  return finalMatrix.rearrangeMatrix(arrangement);
}

void findNodes(Graph & g, set<Node> & s, set<Node> & result)
{
	auto mst=g.minimal_spanning_forest();
	for (const auto& firstNode : s) {
		Neighborhood nextToFirstNode = g.get_neighbors(firstNode).first;
		for (const auto& secondNode : nextToFirstNode) {
			Edge e;
			e.first = firstNode;
			e.second = secondNode;
			if (e.first< e.second && mst.find(e) == mst.end())
			{
				if (s.find(firstNode) != s.end())
				{
					result.insert(firstNode);
				}
				else
				{
					result.insert(secondNode);
				}
			}
		}
	}
}

set<Node> solveDegree3(Graph& g, set<Node>& s, int seed)
{
	set<Node> feedBackSet;
	auto result = graphToMatrix(g, s);
	result.first.print("IncidenceMatrix");
	mat instance = colinearToLinear(result.first);
	instance.print("transformed");
	if (instance==eye<mat>(0,0))
	{
		return feedBackSet;
	}
	Galois ga;
	ga.set_w(16);
	//ga.set_mode_naive();
	ga.set_mode_logtb();
	ga.seed(seed);
	int length;
	int* res = simple_parity_fast(ga, instance.toNMatrix(), instance.getHeight(), instance.getWidth(), &length);
	for (int i = 0; i < length; i++)
	{
		cout << "Delete edge from " << result.second[res[i]].first << " to " << result.second[res[i]].second << endl;
		g.remove_edge(result.second[res[i]].first, result.second[res[i]].second);
	}
	for (int i=0;i<length;i+=2)
	{
		feedBackSet.insert(result.second[res[i]].first);
	}
	findNodes(g, s, feedBackSet);
	for (Node n : feedBackSet)
	{
		cout << "Delete node: " << n << endl;
	}
	delete[] res;
	return feedBackSet;
}


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

