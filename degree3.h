#pragma once
#ifndef debug
    #define debug if(1)
#endif
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include "graph.hpp"
#include "lin_parity.h"
#include "matrix.h"
#include "galois.h"

using namespace FvsGraph;

void generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, mat & matrix, std::size_t columnNumber, const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow);
pair<mat,vector<Edge>> graphToMatrix(const Graph& g, const set<Node>& u, const vector<int>& nodeToComponent);
mat transformFullRowRank(mat input);
mat matrixToStadardForm(const mat & input);
set<Node> solveDegree3(Graph& g, set<Node>& s, const vector<int>& nodeToComponent);
int getComponentNumber(const Node & n,const vector<int>& nodeToComponent);
void findNodes(Graph & g, set<Node> & s, set<Node> & result);
