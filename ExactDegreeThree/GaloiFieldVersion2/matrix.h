#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <functional>
#include "galois.h"

class mat
{
public:
	mat(Galois& ga, int width, int height);
	mat(int width, int height);
	mat();
	mat(uint64_t** pmatrix,int width, int height);
	mat(std::initializer_list<std::initializer_list<uint64_t>> lst);
	~mat();
	void zeros();
	mat operator+(const mat& rhs);
	mat operator*(const mat& rhs);
	bool operator==(const mat & rhs);
	bool operator!=(const mat & rhs);
	uint64_t & operator()(int a, int b);
	uint64_t operator()(int a, int b) const;
	uint64_t at(const int a, const int b) const;
	uint64_t & at(const int a, const int b);
	const int getHeight() const
	{
		if (matrix.empty())
		{
			return 0;
		}
		else
		{
			return matrix[0].size();
		}
	};
	const int getWidth() const 
	{
		return matrix.size();
	};
	void print(const std::string & str) const;
	std::vector<uint64_t> & col(const int colnumber);
	std::vector<uint64_t> col(const int colnumber) const;
	//std::vector<uint64_t> getColumn(const int colNumber) const;
	void addColumn(const std::vector<uint64_t> column);
	void shed_row(const int rowNumber);
	//void deleteRow(const int row);
	int n_rows;
	int n_cols;
	bool is_square() const;
	mat t();
	mat i();
	mat operator-(const mat & rhs);
	std::tuple<mat, mat,std::vector<int>> upper_triangle_transform();
	mat extractColumns(const std::vector<int> & index);
	uint64_t det();
	void swapColumns(int a, int b);
	mat backSubstituation(std::tuple<mat, mat, std::vector<int>>& upperTriangle);
	void extractMatrix(std::vector<int> rowCols);
	std::vector<int> maxSubmatrix();
	uint64_t** toNMatrix();
	mat rearrangeMatrix(const std::vector<int> & arrangement);
//	std::vector<int> fullRankMatrixPosition();
//	std::pair<std::vector<int>, mat> inverseSubmatrix();

private:
	void updateDimension();
	Galois g;
	std::vector<std::vector<uint64_t>> matrix;
	int w;
	int h;
	void columnTransform(mat & matrix, mat & inverse, int colNumber, int destColNumber, int rowNumber);
	void columnOperation(mat& matrix, int colNumber, uint64_t factor);
	void columnOperation(mat & matrix, int firstColNumber, int secondColNumber, uint64_t factor);
	void addValue(uint64_t value,int columnNumber);
	uint64_t** nMat=nullptr;
	void freeMat();
	int nMatHeight = 0;
	//void newRow();
	//int currentRow = 0;
};

/*class column
{
public:
	column(std::vector<uint64_t> & c)
	{
		col = c;
	};
	int getHeight()
	{
		return col.size();
	};
	std::vector<uint64_t> & col;
};*/

mat join_rows(const mat & left,const mat & right);
mat join_rows(const mat & left, const std::vector<uint64_t> & right);
int findNonZero(mat & input, int row, int startCol);
//void lu(mat & l, mat & u, mat & p, mat & input);
int matRank(mat Matrix);
/*template <typename T>
T eye(const int numRow, const int numCol);*/
template <typename T>
T eye(const int numRow, const int numCol)
{
	T id(numRow, numCol);
	for (int i = 0; i<std::min(numRow, numCol); i++)
	{
		id(i, i) = 1;
	}
	return id;
};

std::vector<int> findRedundantRows(uint64_t** pmatrix, int height, int width);
