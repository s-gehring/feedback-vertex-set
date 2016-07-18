#include "degree3.h"

using namespace std;

int getComponentNumber(const Node & n, const vector<int> & nodeToComponent)
{
	return nodeToComponent[n];
}

void generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, mat & matrix, std::size_t columnNumber ,const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow)
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
}

pair<mat,vector<Edge>> graphToMatrix(const Graph& g, const set<Node>& u,const vector<int> & nodeToComponent)
{
  set<pair<Edge,Edge>> edgePairs;
	vector<int> nodeToRow(nodeToComponent.size());
	vector<int> componentToRow(nodeToComponent.size(), -1);
  int lastUsedRow = 0;
  for (const auto &it : g.get_adjacency_list()) {
	  if (u.find(it.first) != u.end())
	  {
		  nodeToRow[it.first] = lastUsedRow;
		  lastUsedRow++;
	  }
	  else
	  {
		  if (componentToRow[getComponentNumber(it.first,nodeToComponent)] == -1)
		  {
			  componentToRow[getComponentNumber(it.first,nodeToComponent)] = lastUsedRow;
			  lastUsedRow++;
		  }
		  nodeToRow[it.first] = componentToRow[getComponentNumber(it.first,nodeToComponent)];
	  }
  }
  cout <<"Number of Nodes: "<<lastUsedRow << endl;
  vector<int> pairNumber(g.get_n(),0);
  size_t edgeNumber = 0;
  for (const auto& firstNode : u) {
    vector<Edge> neighbours;
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
      }
      pairNumber[nodeToRow[firstNode]]+=neighbours.size()-1;
      for (size_t i=0;i<neighbours.size();i++)
      {
        for(size_t j=0;j<i;j++)
        {
            edgePairs.insert(make_pair(neighbours[i],neighbours[j]));
        }
      }
  }
  cout<<"Number of edges: "<<edgeNumber <<" Number of nEdges: "<<g.get_m() <<endl << "Number of pairs: " << edgePairs.size() << endl <<"Max Index: " <<nodeToComponent.size()<< endl;
  vector<Edge> assignment;
  vector<vector<int>> edgesUsed;
  vector<vector<int>> lastVertexIndex;
  for (size_t i=0;i<g.get_n();i++)
  {
    vector<int> v(g.get_n(),0);
    vector<int> w(g.get_n(), -1);
    edgesUsed.push_back(v);
    lastVertexIndex.push_back(w);
  }
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
    assignment.push_back(p.first);
    assignment.push_back(p.second);
  }
  return make_pair(matrix,assignment);
}

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
  mat fullRank = input;
  fullRank.shed_row(0);
  if (fullRank.is_square())
  {
	  return eye<mat>(0, 0);
  }
  auto standardForm = fullRank.toStandarForm();
  vector<int> arrangement = standardForm.second;
  vector<int> restIndex;
  for (size_t i = fullRank.getHeight();i<fullRank.getWidth();i++)
  {
	  restIndex.push_back(i);
  }
  mat newRight = standardForm.first.extractColumns(restIndex);
  mat finalMatrix= join_rows(newRight.t(),eye<mat>(fullRank.n_cols-fullRank.n_rows, fullRank.n_cols-fullRank.n_rows));
  return finalMatrix.rearrangeMatrix(arrangement);
}

void print_edges(const set<Edge>& s) {
    set<Edge>::iterator it = s.begin();
    if (s.size() > 0) {
      cout << "{" << (*it).first<<"-" << (*it).second;
      while (++it != s.end()) {

        cout << ", " << (*it).first<<"-" << (*it).second;
      }
      cout << "}" << endl;
    }
    else {
      cout << "{}" << endl;
    }
  }

void findNodes(Graph & g, set<Node> & s, set<Node> & result)
{
	auto mst=g.minimal_spanning_forest();
    for(const auto &it : g.get_adjacency_list()) {
		Node firstNode=it.first;
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

set<Node> solveDegree3(Graph& g, set<Node>& s, int seed, const vector<int> & nodeToComponent,Galois & ga)
{
  g.print_tidy();
  g.print_nodeset(s);
	set<Node> feedBackSet;
	auto result = graphToMatrix(g, s,nodeToComponent);
	result.first.print("IncidenceMatrix");
	mat instance = colinearToLinear(result.first);
	instance.print("transformed");
	if (instance==eye<mat>(0,0))
	{
		return feedBackSet;
	}
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
