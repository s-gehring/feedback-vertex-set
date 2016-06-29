#include "Tests.h"



Tests::Tests()
{
}


Tests::~Tests()
{
}

void Tests::testAll()
{
	initialize();
	add();
	t();
	multi();
	shed_row();
	i();
	rankMat();
	extractMatrix();
	maxSubmatrix();
	toStandardForm();
	col();
	extractColumns();
	inverseSubmatrix();
}

void Tests::add()
{
	mat m = { {1},{1},{3} };
	mat r(3,1);
	assert(m + m ==r);
}

void Tests::initialize()
{
	mat m = { { 1 },{ 1 },{ 3 } };
	mat r(3, 1);
	r(0, 0) = 1;
	r(1, 0) = 1;
	r(2, 0) = 3;
	assert(m == r);
}

void Tests::t()
{
	mat m = { { 1 },{ 1 },{ 3 } };
	mat n = { { 1,1,3 } };
	assert(m.t() == n);
}

void Tests::multi()
{
	mat m = { {1,0} };
	mat n = { {0},{1} };
	mat r = { {0,0} ,{1,0} };
	assert(n*m == r);
}

void Tests::shed_row()
{
	mat m = { { 0,0 } ,{ 1,0 } };
	mat r = { {0,0}};
	m.shed_row(1);
	assert(m == r);
	assert(m.getHeight()== 1);
}

void Tests::i()
{
	mat m = { {1,1}, {1,0} };
	assert(m*m.i()==eye<mat>(2,2));
	mat m3 = { { 1,1,1 },{ 0,0,1 },{ 0,1,1 } };
	assert(m3*m3.i() == eye<mat>(3, 3));
	mat m2 = { { 23,9,45 },{ 0,0,49 },{0,54,40} };
	assert(m2*m2.i() == eye<mat>(3, 3));
}

void Tests::rankMat()
{
	mat m = { { 1,1 },{ 1,0 } };
	assert(matRank(m) == 2);
	mat m2 = { { 1,1 },{ 2,2 } };
	assert(matRank(m2) == 1);
	mat m3 = { { 1,1,0 },{ 0,0,1 },{1,1,1} };
	assert(matRank(m3) == 2);
	mat m4 = { { 1,1,1,0 },{ 1,0,0,1 },{ 0,1,1,1 } };
	assert(matRank(m4) == 2);
	mat m5 = { { 0,0,1,1,1,0 },{0,0, 1,0,0,1 },{0,0 ,0,1,1,1 } };
	assert(matRank(m4) == 2);
}

void Tests::extractMatrix()
{
	mat m3 = { { 1,1,0 },{ 0,0,1 },{ 1,1,1 } };
	mat r3= { { 1,0 },{ 1,1 } };
	m3.extractMatrix(std::vector<int>({ 1 }), std::vector<int>({0, 1,2 }));
	assert(m3== r3);
	mat m = { { 1,0,1,0 },{1, 0,0,1 },{ 0,1,1,1 },{0,1,2,3} };
	mat r = { { 1,1 },{ 0,1 } };
	m.extractMatrix(std::vector<int>({ 1,3 }), std::vector<int>({ 0,1,2,3 }));
	assert(m == r);
	mat m2 = { { 1,0,1,0 },{ 1, 0,0,1 },{ 0,1,1,1 },{ 0,1,2,3 } };
	mat r2 = { { 0,1,0 },{ 1,1,1 },{1,2,3} };
	m2.extractMatrix(std::vector<int>({ 0}), std::vector<int>({ 1,0,2,3 }));
	assert(m2 == r2);
	mat m4 = { { 1,0,1,0 },{ 1, 0,0,1 },{ 0,1,1,1 },{ 0,1,2,3 } };
	mat r4 = { { 0,1},{ 1,1}};
	m4.extractMatrix(std::vector<int>({ 0,3 }), std::vector<int>({ 1,0,2,3 }));
	assert(m4 == r4);
}

void Tests::maxSubmatrix()
{
	mat m3 = { { 1,0,5 },{ 0,0,0 },{ 1,0,1 } };
	assert(m3.maxSubmatrix() == std::vector<int>({1}));
	mat m = { { 1,0,1,0 },{ 0, 0,0,0 },{ 0,0,1,0 },{ 0,0,0,0 } };
	assert(m.maxSubmatrix() == std::vector<int>({ 3,1 }));
	mat m2 = { { 1,1,1},{ 0,0,1},{ 0,0,1}};
	assert(m2.maxSubmatrix() == std::vector<int>({ 1 }));
}

void Tests::extractColumns()
{
	mat m3 = { { 1,1,0 },{ 0,0,1 },{ 1,1,1 } };
	mat r3 = { { 1,0 },{ 0,1 }, {1,1} };
	assert(m3.extractColumns(std::vector<int>({ 0,2 })) == r3);
	mat m = { { 1,0,1,0 },{ 1, 0,0,1 },{ 0,1,1,1 },{ 0,1,2,3 } };
	mat r = { {1},{1},{0},{0} };
	assert(m.extractColumns(std::vector<int>({ 0 })) == r);
}

void Tests::col()
{
	mat m3 = { { 1,1,0 },{ 0,0,1 },{ 1,1,1 } };
	assert(m3.col(2) == std::vector<uint64_t>({0,1,1 }));
}

std::vector<int> Tests::findRest(const std::vector <int>& index, int dim)
{
	std::vector<int> result;
	for (int i = 0; i<dim; i++)
	{
		if (std::find(index.begin(), index.end(), i) == index.end())
		{
			result.push_back(i);
		}
	}
	return result;
}

void Tests::inverseSubmatrix()
{
	/*mat m3 = { { 1,1,0 },{ 0,0,1 },{ 0,1,0 } };
	assert(m3.inverseSubmatrix().first == std::vector<int>({ 0,2,1 }));
	assert(m3.inverseSubmatrix().second == m3.i());
	mat m = { {0, 1,1,0 },{0, 0,0,1 },{ 0,0,1,0 } };
	assert(m.inverseSubmatrix().first == std::vector<int>({ 1,3, 2}));
	assert(m.inverseSubmatrix().second == m3.i());
	mat m2 = { {1,1,0,0 },{0,0,0,1 },{ 0,1,0,0 } };
	assert(m2.inverseSubmatrix().first == std::vector<int>({ 0,3, 1 }));
	assert(m2.inverseSubmatrix().second == m3.i());
	mat m4 = { { 1,1,0,0 },{ 1,1,0,1 },{ 1,1,1,0 } };
	mat r4 = { { 1,0,0 },{ 1,0,1 },{ 1,1,0 } };
	assert(m4.inverseSubmatrix().second == r4.i());*/

	/*std::vector<int> invSubMatrix = std::get<2>(m4.upper_triangle_transform());
	std::vector<int> restIndex = findRest(invSubMatrix, m4.n_cols);
	invSubMatrix.resize(m4.getHeight());
	mat submatrix = m4.extractColumns(invSubMatrix);
	submatrix.print("sub");
	mat inverse=submatrix.i();
	r4.i().print("r4i");
	assert(inverse == r4.i());*/
}

void Tests::toStandardForm()
{
	mat m4 = { { 1,1,0,0 },{ 1,1,0,1 },{ 0,0,1,0 } };
	mat r4 = { { 1,0,0,1 },{ 0,1,0,0 },{ 0,0,1,0 } };
	//r4.toStandarForm().first.print("r4");
	assert(m4.toStandarForm().first == r4);

	/*mat m3 = { { 1,0,5 },{ 1,0,0 }};
	std::cout << "------" << std::endl;
	for (int a : m3.fullRankMatrixPosition())
	{
		std::cout << a << std::endl;
	}
	assert(m3.fullRankMatrixPosition() == std::vector<int>({ 0,2 }));
	mat m = { { 1,0,1,0 },{ 0, 0,0,0 },{ 0,0,1,0 },{ 0,0,0,0 } };
	assert(m.fullRankMatrixPosition() == std::vector<int>({ 0,1 }));
	mat m2 = { { 1,1,1 },{ 0,0,1 },{ 0,0,1 } };
	assert(m2.fullRankMatrixPosition() == std::vector<int>({ 0,2 }));*/
}
