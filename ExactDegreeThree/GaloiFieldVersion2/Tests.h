#pragma once
#include <cassert>
#include "matrix.h"

class Tests
{
public:
	Tests();
	~Tests();
	void testAll();
	void add();
	void initialize();
	void t();
	void multi();
	void i();
	void shed_row();
	void rankMat();
	void extractMatrix();
	void maxSubmatrix();
	void toStandardForm();
	void extractColumns();
	void col();
	void inverseSubmatrix();
	std::vector<int> findRest(const std::vector <int>& index, int dim);
};

