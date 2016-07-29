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
	mat(int height, int width);
	mat();
	mat(uint64_t** pmatrix,int height, int width);
	mat(std::initializer_list<std::initializer_list<uint64_t>> lst);
	void zeros();
	mat operator+(const mat& rhs);
	mat operator*(const mat& rhs);
	bool operator==(const mat & rhs);
	bool operator!=(const mat & rhs);
	uint64_t & operator()(const int a,const int b);
	uint64_t operator()(const int a,const int b) const;
	uint64_t at(const int a, const int b) const;
	uint64_t & at(const int a, const int b);
	std::size_t getHeight() const
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
	std::size_t getWidth() const 
	{
		return matrix.size();
	};
	void print(const std::string & str) const;
	std::vector<uint64_t> & col(const int colnumber);
	std::vector<uint64_t> col(const int colnumber) const;
	void addColumn(const std::vector<uint64_t> & column);
	void shed_row(const int rowNumber);
	std::size_t n_rows;
	std::size_t n_cols;
	bool is_square() const;
	mat t();
	mat i();
	mat operator-(const mat & rhs);
	std::tuple<mat, mat,std::vector<int>> upper_triangle_transform();
	mat extractColumns(const std::vector<int> & index);
	uint64_t det();
	void swapColumns(int a, int b);
	mat backSubstituation(std::tuple<mat, mat, std::vector<int>>& upperTriangle);
	std::vector<int> maxSubmatrix();
	uint64_t** toNMatrix();
	mat rearrangeMatrix(const std::vector<int> & arrangement);
	std::pair<mat, std::vector<int>> toStandarForm();

private:
	void updateDimension();
	std::vector<std::vector<uint64_t>> matrix;
	void rowTransform(mat & matrix, mat & inverse, int rowNumber, int destRowNumber, int colNumber);
	void rowOperation(mat& matrix, int rowNumber, uint64_t factor);
	void rowOperation(mat & matrix, int firstRowNumber, int secondRowNumber, uint64_t factor);
	void rowTransform(mat & matrix, int rowNumber, int destRowNumber, int colNumber);

	void columnTransform(mat & matrix, mat & inverse, int colNumber, int destColNumber, int rowNumber);
	void columnOperation(mat& matrix, int colNumber, uint64_t factor);
	void columnOperation(mat & matrix, int firstColNumber, int secondColNumber, uint64_t factor);
	void addValue(uint64_t value,int columnNumber);
};

mat join_rows(const mat & left,const mat & right);
mat join_rows(const mat & left, const std::vector<uint64_t> & right);
int findNonZero(mat & input, int row, int startCol);
int matRank(mat Matrix);
template <typename T>
T eye(const int numRow, const int numCol)
{
	T id(numRow, numCol);
	for (int i = 0; i<std::min(numRow, numCol); i++)
	{
		id(i, i) = 1;
	}
	return id;
}
std::vector<int> findRedundantRows(uint64_t** pmatrix, int height, int width);