#include "degree3.h"

using namespace std;

 /**
	* @brief generates the incidence of an edge
	* @param [in] g input graph
  * @param [in] e set of nodes can be included in the final fvs set
	* @param [in] edgesUsed the number of artificial edges we already inserted for a specific edge
	* @param [in] lastUsedRow the last index of the row we added to the matrix
	* @param [in] matrix the resulting matrix
	* @param [in] columnNumber the number of columns we already generated in matrix
	* @param [in] pairNumber the number of edge pairs that shares a specific node
	* @param [in] lastVertexIndex the row of the last node of artificial edges for a specific edge, -1 if none was generated so far
	* @param [in] nodeToRow assigns a node to its initial row number
	*/
void generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, mat & matrix, std::size_t columnNumber ,const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow)
{
	int mi = min(nodeToRow[e.first], nodeToRow[e.second]);
	int ma = max(nodeToRow[e.first], nodeToRow[e.second]);
	edgesUsed[mi][ma]++;
	//set the entry of the first node to 1
	if (lastVertexIndex[mi][ma] != -1)
	{
		matrix(lastVertexIndex[mi][ma], columnNumber) = 1;
	}
	else
	{
		matrix(mi, columnNumber) = 1;
	}
	//set the entry of the second node to 1
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

 /**
	* @brief generates a colinear matroid parity problem from the graph
	* @param [in] g input graph
  * @param [in] s set of nodes can be included in the final fvs set
	* @param [in] nodeToComponent assigns all nodes not in s to an connected component number
	* @returns the matrix and the assignment of column to the corresponding edge
	*/
pair<mat,vector<Edge>> graphToMatrix(const Graph& g, const set<Node>& s,const vector<int> & nodeToComponent)
{
	//first find all edge pairs that shares a node in s
  set<pair<Edge,Edge>> edgePairs;
	vector<int> nodeToRow(nodeToComponent.size());
	vector<int> componentToRow(nodeToComponent.size(), -1);
  int lastUsedRow = 0;
  for (const auto &it : g.get_adjacency_list()) {
	  if (s.find(it.first) != s.end())
	  {
		  nodeToRow[it.first] = lastUsedRow;
		  lastUsedRow++;
	  }
	  else
	  {
		  if (componentToRow[nodeToComponent[it.first]] == -1)
		  {
			  componentToRow[nodeToComponent[it.first]] = lastUsedRow;
			  lastUsedRow++;
		  }
		  nodeToRow[it.first] = componentToRow[nodeToComponent[it.first]];
	  }
  }
  vector<int> pairNumber(g.get_n(),0);
  size_t edgeNumber = 0;
  for (const auto& firstNode : s) {
    vector<Edge> neighbours;
	  Neighborhood nextToFirstNode = g.get_neighbors(firstNode).first;
	  for (const auto& secondNode : nextToFirstNode) {
		  if (s.find(secondNode)==s.end() || firstNode < secondNode)
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
	//create the result matrix
  int row_number = lastUsedRow + 2 * edgePairs.size() - edgeNumber;
  mat matrix(row_number,edgePairs.size()*2);
  int columnNumber=0;
	//generate the incidence vector for the edges
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
 /**
	* @brief converts the colinear matroid parity problem to a linear parity problem
	*/
mat colinearToLinear(const mat & input)
{
	//delete one row to get a full rank matrix since graph is connected so the matrix has rank row+1
  mat fullRank = input;
  fullRank.shed_row(0);
	//if the matrixx is a square, the graph is a tree
  if (fullRank.is_square())
  {
	  return eye<mat>(0, 0);
  }
	//transform the matrix to the from (IB) where I is the identity matrix
  auto standardForm = fullRank.toStandarForm();
  vector<int> arrangement = standardForm.second;
  vector<int> restIndex;
  for (size_t i = fullRank.getHeight();i<fullRank.getWidth();i++)
  {
	  restIndex.push_back(i);
  }
	//extract the right part B of the standard form matrix
  mat newRight = standardForm.first.extractColumns(restIndex);
	//transpose and concat with identity matrix
  mat finalMatrix= join_rows(newRight.t(),eye<mat>(fullRank.n_cols-fullRank.n_rows, fullRank.n_cols-fullRank.n_rows));
  return finalMatrix.rearrangeMatrix(arrangement);
}

 /**
	* @brief find remaining nodes of the fvs, that is all the nodes which are not a common node of a two group 
	* @param [in] g input graph
  * @param [in] s set of nodes can be included in the final fvs set
	* @param [in] fvs set
	*/
void findNodes(Graph & g, set<Node> & s, set<Node> & result)
{
		auto mst=g.minimal_spanning_forest();
		//deletes all nodes that are not in the mst
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
 /**
	* @brief solves the fvs problem in the degree 3 case in polynomial time
	* @param [in] g input graph it should be connected without vertices of degree 1
  * @param [in] s set of nodes can be included in the final fvs set
	* @param [in] nodeToComponent an array that assigns all nodes not in s to an connected component number
	* @returns fvs set
	*/
set<Node> solveDegree3(Graph& g, set<Node>& s, const vector<int> & nodeToComponent)
{
	set<Node> feedBackSet;
	//convert graph to matrix
	auto result = graphToMatrix(g, s,nodeToComponent);
	//convert colinear parity problem to linear parity problem
	mat instance = colinearToLinear(result.first);
	//return emtpy fvs if graph is a tree
	if (instance==eye<mat>(0,0))
	{
		return feedBackSet;
	}
	int length;
	//solve the linear parity priblem
	int* res = simple_parity_fast(Galois::getInstance(), instance.toNMatrix(), instance.getHeight(), instance.getWidth(), &length);
	//remove the edges of the solution of the colinear parity problem
	for (int i = 0; i < length; i++)
	{
		g.remove_edge(result.second[res[i]].first, result.second[res[i]].second);
	}
	//insert common node of the edgepairs to the fvs
	for (int i=0;i<length;i+=2)
	{
		feedBackSet.insert(result.second[res[i]].first);
	}
	//find the remaining nodes of the fvs
	findNodes(g, s, feedBackSet);
	delete[] res;
	return feedBackSet;
}
