#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <ctime>
#include "Tests.h"
#include "lin_parity.h"
#include "matrix.h"
#include "graph.hpp"
#include "fvs_solver.hpp"


//vector<uint64_t> generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, int row_number, const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow);
void generateIncidenceVector(const Graph& g, const Edge & e, vector<vector<int>>& edgesUsed, int & lastUsedRow, mat & matrix, int columnNumber, const vector<int> & pairNumber, vector<vector<int>>& lastVertexIndex, vector<int> & nodeToRow);
pair<mat,vector<Edge>> graphToMatrix(const Graph& g, const set<Node>& u);
//vector<int> findRest(const vector <int>& index, int dim);
//mat rearrangeMatrix(const mat & input, const vector<int> & arrangement);
mat transformFullRowRank(mat input);
mat matrixToStadardForm(const mat & input);
set<Node> solveDegree3(Graph& g, set<Node>& s,int seed);
void findNodes(Graph & g, set<Node> & s, set<Node> & result);
