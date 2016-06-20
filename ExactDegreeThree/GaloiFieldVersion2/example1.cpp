#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <vector>
#include <algorithm>
#include <tuple>
#include <ctime>
#include "fvs_solver.h"
#include "Tests.h"
//#include "lin_parity.h"
#include "matrix.h"

using namespace std;
using namespace fvs;

mat generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow,int row_number, const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex)
{
  mat result(row_number,1);
  result.zeros();
  int mi=min(source(e,g),target(e,g));
  int ma=max(source(e,g),target(e,g));
  edgesUsed[mi][ma]++;
  result(lastVertexIndex[mi][ma],0)=1;
  //cout<<mi<<" "<<ma<<" "<<lastVertexIndex[mi][ma]<<endl;
  if (edgesUsed[mi][ma]==pairNumber[mi]+pairNumber[ma])
  {
    //result(ma,0)=-1;
	  result(ma, 0) = 1;
	//cout<<"case 1 "<<ma<<endl;

  }
  else
  {
    //result(lastUsedRow,0)=-1;
	result(lastUsedRow, 0)= 1;
    lastVertexIndex[mi][ma]=lastUsedRow;
    lastUsedRow++;
  }
  return result;
}

pair<mat,vector<Edge>> graphToMatrix(const Graph& g, const set<Node>& u)
{
  typedef graph_traits<Graph>::vertices_size_type id;

  typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
  set<pair<Edge,Edge>> edgePairs;
  typedef graph_traits<Graph>::vertex_iterator node_iterator;

  map<Node,int> nodeIndex;
  int index=0;
  
  /*pair<node_iterator, node_iterator> nIt = vertices(g);
  for (node_iterator it = nIt.first; it != nIt.second; ++it)
  {
    nodeIndex[*it]=index;
    index++;
  }*/
  vector<int> pairNumber(num_vertices(g),0);
  for (const auto& n : u) {
      size_t edgesToOtherUs = 0;
      vector<Edge> neighbours;
      pair<edge_iterator, edge_iterator> eIt = out_edges(n, g);
      for (edge_iterator it = eIt.first; it != eIt.second; ++it) {
        neighbours.push_back(*it);
        //edgePairs.insert(make_pair(nodeIndex[i],nodeIndex[target(*it, g)]));
        //neighbours.push_back(nodeIndex[target(*it, g)]);
      }
      pairNumber[n]=neighbours.size()-1;
      for (int i=0;i<neighbours.size();i++)
      {
        for(int j=0;j<i;j++)
        {
            edgePairs.insert(make_pair(neighbours[i],neighbours[j]));
        }
      }
  }

  mat matrix;
  vector<Edge> assignment;
  int row_number=num_vertices(g)+2*edgePairs.size()-num_edges(g);
  vector<vector<int>> edgesUsed;
  vector<vector<int>> lastVertexIndex;
  for (int i=0;i<num_vertices(g);i++)
  {
    vector<int> v(num_vertices(g),0);
    vector<int> w(num_vertices(g),i);
    edgesUsed.push_back(v);
    lastVertexIndex.push_back(w);
  }
  int lastUsedRow=num_vertices(g);
  /*for (int i=0;i<num_vertices(g);i++)
  {
    lastVertexIndex[i]=i;
  }*/
  for(auto & p: edgePairs)
  {
    mat first=generateIncidenceVector(g,p.first,edgesUsed,lastUsedRow,row_number,pairNumber,lastVertexIndex);
    mat second=generateIncidenceVector(g,p.second,edgesUsed,lastUsedRow,row_number,pairNumber,lastVertexIndex);
    matrix=join_rows(matrix,first);
    matrix=join_rows(matrix,second);
    assignment.push_back(p.first);
    assignment.push_back(p.second);
  }
  return make_pair(matrix,assignment);
}

vector<int> findRest(const vector <int>& index, int dim)
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
}

mat extractColumns(const mat & input, const vector<int> & index)
{
    mat result;
    for(int i: index)
    {
      result=join_rows(result,input.col(i));
    }
    return result;
}

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

mat rearrangeMatrix(const mat & input, const vector<int> & arrangement)
{
  mat result(input.n_rows,input.n_cols);
  for(int i=0;i<arrangement.size();i++)
  {
    result.col(arrangement[i])=input.col(i);
  }
  return result;
}


mat transformFullRowRank(mat input)
{
  int i=0;
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
  }
  return input;
}

mat matrixToStadardForm(const mat & input)
{
  mat fullRank= transformFullRowRank(input);
  fullRank.print("transformed row rank");
  if (fullRank.is_square())
  {
    return eye<mat>(fullRank.n_rows,fullRank.n_cols);
  }
  /*mat L, U, P;
  lu(L,U,P,fullRank.t());
  L.print("L");
  U.print("U");
  P.print("P");
  mat res=U.t()*L.t()*P;
  res.print("check");
  vector<int> invSubMatrix= findOnes(P,L.n_cols);*/
  vector<int> invSubMatrix =std::get<2>(fullRank.upper_triangle_transform());
  invSubMatrix.resize(fullRank.getHeight());
  vector<int> restIndex= findRest(invSubMatrix,fullRank.n_cols);
  mat submatrix=extractColumns(fullRank,invSubMatrix);
  submatrix.print("sub");
  mat inverse=submatrix.i();
  inverse.print("inverse");
  mat rest= extractColumns(fullRank,restIndex);
  rest.print("rest");
  mat newRight=inverse*rest;
  newRight.print("newRight");
  vector<int> arrangement=invSubMatrix;
  arrangement.insert(arrangement.end(), restIndex.begin(), restIndex.end() );
  mat finalMatrix= join_rows(newRight.t(),eye<mat>(fullRank.n_cols-fullRank.n_rows, fullRank.n_cols-fullRank.n_rows));
  finalMatrix.print("finalMatrix");
  return rearrangeMatrix(finalMatrix,arrangement);
}

int main(int argc, char** argv)
  {
  std::srand(std::time(0));
  Graph g;
  read_graph(g, "mini_graph.txt");
  print_graph(g);
  typedef graph_traits<Graph>::vertex_iterator node_iterator;
  pair<node_iterator, node_iterator> nIt = vertices(g);
  set<Node> s;
  for (node_iterator it = nIt.first; it != nIt.second; ++it)
  {
    s.insert(*it);
  }
  Tests test;
  test.testAll();
  auto result=graphToMatrix(g,s);
  result.first.print("IncidenceMatrix");
  mat instance= matrixToStadardForm(result.first);
  instance.print("transformed");
  /*vector<int> test = simple_parity(instance, 1000);
  for(int i: test)
  {
      cout<<"Delete edge from "<<source(result.second[i],g)<< " to "<<target(result.second[i],g) <<endl;
  }*/
   //cout << "Final parity basis are the columns:"<<endl << test << endl;
  /*for (int i=0;i<vec.n_elem;i++)
  {
    int col1=(int) (vec(i)+0.2);
    int col2=(int) (vec(i)*0.2); 
    cout<<"Delete edge from "<<source(result.second[col1],g)<< " to "<<target(result.second[col2],g) <<endl;
  };*/
  /*for(int i=0;i<test.size();i+=2)
  {
      cout<<"Delete note from "<<source(result.second[i],g) <<endl;
      remove_vertex(source(result.second[i],g),g);
  }
  print_graph(g);*/
  return 0;
  }

