#include "matrix.h"
#include "galois.h"


mat::mat(Galois& ga,int height, int width)
{
	g=ga;
	for (int i = 0; i<width; i++)
	{
		std::vector<uint64_t>column(height, 0);
		matrix.push_back(column);
	}
	updateDimension();
}

mat::mat(uint64_t** pmatrix, int height, int width) 
{
	g.set_w(64);
	g.set_mode_naive();
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			addValue(pmatrix[i][j],j);
		}
		//newRow();
	}
}


mat::mat(int height, int width)
{
	g.set_w(64);
	g.set_mode_naive();
	for (int i = 0; i<width; i++)
	{
		std::vector<uint64_t>column(height,0);
		matrix.push_back(column);
	}
	updateDimension();
}

mat::mat()
{
	g.set_w(64);
	g.set_mode_naive();
	updateDimension();
}

mat::~mat()
{
	freeMat();
}

mat::mat(std::initializer_list<std::initializer_list<uint64_t>> lst)
	: n_rows(lst.size()), n_cols(n_rows ? lst.begin()->size() : 0)
{
	g.set_w(64);
	g.set_mode_naive();
	int i = 0;
	for (const auto& row : lst) {
		for (const auto& value : row) {
			addValue(value,i);
			i++;
		}
		i=0;
		//newRow();
	}
}

void mat::addValue(uint64_t value, int columnNumber)
{
	if (matrix.size()<= columnNumber)
	{
		matrix.push_back(std::vector<uint64_t>(1, value));
	}
	else
	{
		matrix[columnNumber].push_back(value);
	}
}

/*void mat::newRow()
{
	currentRow = 0;
}*/

void mat::zeros()
{
	for (auto & i: matrix)
	{
		for (auto & j : i)
		{
			j = 0;
		}
	}
}

uint64_t & mat::operator()(const int a,const int b)
{
	return matrix[b][a];
}

uint64_t mat::operator()(const int a, const int b) const
{
	return matrix[b][a];
}

uint64_t mat::at(const int a, const int b) const
{
	return matrix[b][a];
}

uint64_t & mat::at(const int a, const int b)
{
	return matrix[b][a];
}

mat mat::operator*(const mat & rhs)
{
	if (getWidth() != rhs.getHeight())
	{
		throw;
	}
	mat result(getHeight(), rhs.getWidth());
	for (int i = 0; i < result.getWidth(); i++)
	{
		for (int j=0;j<result.getHeight();j++)
		{
			for (int k=0;k<getWidth();k++)
			{
				result(j,i) = g.add(result(j,i),g.multiply(matrix[k][j],rhs(k,i)));
			}
		}
	}
	return result;
}

mat mat::operator+(const mat & rhs)
{
	if (getWidth() != rhs.getWidth() || getHeight()!= rhs.getHeight())
	{
		throw;
	}
	mat result(getHeight(), getWidth());
	for (int i = 0; i < result.getWidth(); i++)
	{
		for (int j = 0; j<result.getHeight(); j++)
		{
			result(j,i) = g.add(matrix[i][j], rhs(j,i));
		}
	}
	return result;
}

mat mat::operator-(const mat & rhs)
{
	return (*this)+rhs;
}

mat join_rows(const mat & left, const mat & right)
{
	if (left.getHeight() != 0 && right.getHeight() != 0 && left.getHeight() != right.getHeight())
	{
		throw;
	}
	mat result=left;
	for (int i=0;i<right.getWidth();i++)
	{
		result.addColumn(right.col(i));
	}
	return result;
}

mat join_rows(const mat & left, const std::vector<uint64_t> & right)
{
	if (left.getHeight() != 0 && right.size() != 0 && left.getHeight() != right.size())
	{
		throw;
	}
	mat result = left;
	result.addColumn(right);
	return result;
}

void join_rows_fast(mat & left, const std::vector<uint64_t> & right)
{
	if (left.getHeight() != 0 && right.size() != 0 && left.getHeight() != right.size())
	{
		throw;
	}
	left.addColumn(right);
}

void mat::updateDimension()
{
	n_rows = getHeight();
	n_cols = getWidth();
}

std::vector<uint64_t> & mat::col(const int colnumber)
{
	return matrix[colnumber];
}

std::vector<uint64_t> mat::col(const int colnumber) const
{
	return matrix[colnumber];
}

/*mat mat::shed_row(const int row)
{
	mat result = *this;
	result.deleteRow(row);
	return result;
}

void mat::deleteRow(int rowNumber)
{
	for (auto & row: matrix)
	{
		row.erase(row.begin() + rowNumber);
		row.shrink_to_fit();
	}
	updateDimension();
}*/
void mat::shed_row(const int rowNumber)
{
	for (auto & row : matrix)
	{
		row.erase(row.begin() + rowNumber);
		row.shrink_to_fit();
	}
	updateDimension();
}


void mat::addColumn(const std::vector<uint64_t> & column)
{
	matrix.push_back(column);
	updateDimension();
}

/*std::vector<uint64_t> mat::getColumn(const int colNumber) const
{
	if (colNumber >= matrix.size())
	{
		throw;
	}
	return matrix[colNumber];
}*/


void mat::print(const std::string & str) const
{
	std::cout << str<<std::endl;
	if (getWidth()>30 || getHeight()>30)
	{
		std::cout<<"Matrix too big, Width: "<<getWidth()<<" Height: "<<getHeight() <<std::endl;
	}
	else
	{
		std::stringstream output;
		for (int j = 0; j<getHeight(); j++)
		{
			output.str(std::string());
			for (int i = 0; i < getWidth(); i++)
			{
				output<<matrix[i][j] <<" ";
			}
			std::cout<<output.str() << std::endl;
		}
	}
}

int matRank(mat matrix)
{
	auto transform = matrix.upper_triangle_transform();
	int rank = 0;
	for (int i = 0; i < std::min(std::get<0>(transform).getWidth(),std::get<0>(transform).getHeight()); i++)
	{
		if (std::get<0>(transform)(i,i)!=0)
		{
			rank++;
		}
	}
	return rank;
}

bool mat::is_square() const
{
	return getWidth() == getHeight();
}

mat mat::t()
{
	mat result(getWidth(), getHeight());
	for (int i = 0;i<getWidth();i++)
	{
		for (int j=0;j<getHeight();j++)
		{
			result(i,j) = matrix[i][j];
		}
	}
	return result;
}

mat mat::i()
{
	if (!is_square())
	{
		throw;
	}
	auto uTriangle = upper_triangle_transform();
	return backSubstituation(uTriangle);
}

int findNonZero(mat & input,int row, int startCol)
{
	for (int i = startCol;i<input.getWidth();i++)
	{
		if (input(row, i) != 0)
		{
			return i;
		}
	}
	return -1;
}

/*void lu(mat & l, mat & u, mat & p, mat & input)
{
	l = input;
	int dim = std::max(l.getWidth(), l.getHeight());
	p = eye<mat>(dim, dim);
	for (int n = 0; n < std::min(l.getHeight(), l.getWidth()); n++)
	{
		int k = findNonZero(l, n, n);
		if (k!=-1) //or column is 0
		{
			if (k != n)
			{
				l.swapColumns(k, n);
				p.swapColumns(k, n);
			}
		}
	}
}*/

/*std::vector<int> mat::fullRankMatrixPosition()
{
	if (getHeight() > getWidth())
	{
		throw;
	}
	mat A = (*this);
	//int dim = std::max(A.getWidth(), A.getHeight());
	//p = eye<mat>(dim, dim);
	std::vector<int> result;
	for (int i = 0; i < A.getWidth();i++)
	{
		result.push_back(i);
	}
	for (int n = 0; n < std::min(A.getHeight(), A.getWidth()); n++)
	{
		int k = findNonZero(A, n, n);
		if (k != -1) //or column is 0
		{
			if (k != n)
			{
				A.swapColumns(k, n);
				std::swap(result[k], result[n]);
			}
		}
	}
	result.resize(A.getHeight());
	return result;
}*/

std::vector<int> mat::maxSubmatrix()
{
	if (getWidth() != getHeight())
	{
		throw;
	}
	auto uTriangle = upper_triangle_transform();
	std::vector<int> result;
	for (int i = std::min(std::get<0>(uTriangle).getHeight(), std::get<0>(uTriangle).getWidth())-1; i >=0; i--)
	{
		if (std::get<0>(uTriangle)(i, i) == 0)
		{
			result.push_back(std::get<2>(uTriangle)[i]);
		}
	}
	return result;
}

void mat::extractMatrix(std::vector<int> rowCols, const std::vector<int> & arrangement)
{
	std::sort(rowCols.begin(), rowCols.end(), std::greater<int>());
	for (auto i : rowCols)
	{
		matrix.erase(matrix.begin()+i);
	}
	std::sort(rowCols.begin(), rowCols.end(), [&arrangement](int & left, int & right) {
		return arrangement[left] > arrangement[right];
	});
	for (auto i: rowCols)
	{
		for (auto & j : matrix)
		{
			j.erase(j.begin()+ arrangement[i]);
		}
	}
	matrix.shrink_to_fit();
	matrix[0].shrink_to_fit();
	updateDimension();
}

/*template <typename T>
T eye(const int numRow, const int numCol)
{
	T id(numRow, numCol);
	for (int i=0;i<min(numRow,numCol);i++)
	{
		id(i,i) = 1;
	}
	return id;
}*/

void mat::columnTransform(mat & matrix, mat & inverse, int colNumber, int destColNumber, int rowNumber)
{
	if (matrix(rowNumber, colNumber) != 1)
	{
		uint64_t inv=g.inverse(matrix(rowNumber, colNumber));
		columnOperation(matrix,colNumber, inv);
		columnOperation(inverse,colNumber, inv);
	}
	uint64_t fac= matrix(rowNumber,destColNumber);
	columnOperation(matrix, colNumber, destColNumber, fac);
	columnOperation(inverse,colNumber, destColNumber, fac);
}

void mat::columnOperation(mat& matrix, int colNumber, uint64_t factor)
{
	for (int i = 0; i < matrix.getHeight(); i++)
	{
		matrix(i,colNumber) = g.multiply(matrix(i,colNumber),factor);
	}
}

void mat::columnOperation(mat & matrix, int firstColNumber, int secondColNumber, uint64_t factor)
{
	for (int i = 0; i < matrix.getHeight(); i++)
	{
		matrix(i,secondColNumber) = g.add(g.multiply(matrix(i,firstColNumber), factor), matrix(i,secondColNumber));
	}
}

void mat::rowTransform(mat & matrix, mat & inverse, int rowNumber, int destRowNumber, int colNumber)
{
	if (matrix(rowNumber, colNumber) != 1)
	{
		uint64_t inv = g.inverse(matrix(rowNumber, colNumber));
		rowOperation(matrix, rowNumber, inv);
		rowOperation(inverse, rowNumber, inv);
	}
	uint64_t fac = matrix(destRowNumber, colNumber);
	rowOperation(matrix, rowNumber, destRowNumber, fac);
	rowOperation(inverse, rowNumber, destRowNumber, fac);
}

void mat::rowTransform(mat & matrix, int rowNumber, int destRowNumber, int colNumber)
{
	if (matrix(rowNumber, colNumber) != 1)
	{
		uint64_t inv = g.inverse(matrix(rowNumber, colNumber));
		rowOperation(matrix, rowNumber, inv);
	}
	uint64_t fac = matrix(destRowNumber, colNumber);
	rowOperation(matrix, rowNumber, destRowNumber, fac);
}

void mat::rowOperation(mat& matrix, int rowNumber, uint64_t factor)
{
	for (int i = 0; i < matrix.getWidth(); i++)
	{
		if (factor == 0)
		{
			matrix(rowNumber, i) = 0;
		}
		else if (factor!=1)
		{
			matrix(rowNumber, i) = g.multiply(matrix(rowNumber, i), factor);
		}
	}
}

void mat::rowOperation(mat & matrix, int firstRowNumber, int secondRowNumber, uint64_t factor)
{
	for (int i = 0; i < matrix.getWidth(); i++)
	{
		if (factor == 1)
		{
			matrix(secondRowNumber, i) = g.add(matrix(firstRowNumber, i), matrix(secondRowNumber, i));
		}
		else if (factor != 0)
		{
			matrix(secondRowNumber, i) = g.add(g.multiply(matrix(firstRowNumber, i), factor), matrix(secondRowNumber, i));
		}
	}
}

uint64_t mat::det()
{
	auto transform = upper_triangle_transform();
	uint64_t result = std::get<0>(transform)(0,0);
	for (int i = 1; i < std::get<0>(transform).getWidth(); i++)
	{
		result = g.multiply(result, std::get<0>(transform)(i,i));
	}
	return result;
}

void mat::swapColumns(int a, int b)
{
	auto temp = matrix[a];
	matrix[a] = matrix[b];
	matrix[b]=temp;
}

mat mat::backSubstituation(std::tuple < mat, mat,std::vector<int> > & upperTriangle)
{
	for (int i=std::get<0>(upperTriangle).getHeight()-1;i>=0;i--)
	{
		for (int j=0;j<i;j++)
		{
			columnTransform(std::get<0>(upperTriangle), std::get<1>(upperTriangle), i, j, i);
		}
	}
	assert(std::get<0>(upperTriangle) == eye<mat>(std::get<0>(upperTriangle).getHeight(), std::get<0>(upperTriangle).getWidth()));
	return std::get<1>(upperTriangle);
}

std::tuple<mat,mat,std::vector<int>> mat::upper_triangle_transform()
{
	mat matrix = (*this);
	//mat pInverse = eye<mat>(getHeight(), getWidth());
	mat pInverse = eye<mat>(getWidth(), getWidth());
	std::vector<int> result;
	for (int i = 0; i < getWidth(); i++)
	{
		result.push_back(i);
	}
	for (int n = 0; n < std::min(matrix.getHeight(),matrix.getWidth()); n++)
	{
		int nonzero = findNonZero(matrix, n, n);

		if (nonzero!=-1) //or column is 0
		{
			if (nonzero != n)
			{
				matrix.swapColumns(nonzero, n);
				pInverse.swapColumns(nonzero, n);
				std::swap(result[nonzero],result[n]);
			}

			for (int i = n + 1; i < matrix.getWidth(); i++)
			{
				columnTransform(matrix, pInverse, n, i, n);
			}
		}
	}
	return std::make_tuple(matrix,pInverse,result);
}

std::pair<mat,std::vector<int>> mat::toStandarForm()
{
	mat matrix = (*this);
	std::vector<int> result;
	for (int i = 0; i < getWidth(); i++)
	{
		result.push_back(i);
	}
	for (int n = 0; n < std::min(matrix.getHeight(), matrix.getWidth()); n++)
	{
		int nonzero = findNonZero(matrix, n, n);

		if (nonzero != -1) //or column is 0
		{
			if (nonzero != n)
			{
				matrix.swapColumns(nonzero, n);
				std::swap(result[nonzero], result[n]);
			}

			for (int i = n + 1; i < matrix.getHeight(); i++)
			{
				rowTransform(matrix, n, i, n);
			}
		}
	}
	for (int i = matrix.getHeight() - 1; i >= 0; i--)
	{
		for (int j = 0; j<i; j++)
		{
			rowTransform(matrix, i, j, i);
		}
	}
	return std::make_pair(matrix, result);
}

mat mat::extractColumns(const std::vector<int> & index)
{
	mat result;
	for (int i : index)
	{
		result = join_rows(result, col(i));
	}
	return result;
}

bool mat::operator==(const mat & rhs) 
{
	if (getHeight() != rhs.getHeight() || getWidth() != rhs.getWidth())
	{
		return false;
	}
	for (int i=0;i<getHeight();i++)
	{
		for (int j=0;j<getWidth();j++)
		{
			if (matrix[j][i] != rhs(i,j))
			{
				return false;
			}
		}
	}
	return true;
}

bool mat::operator!=(const mat & rhs)
{
	return !((*this) == rhs);
}

std::vector<int> findRedundantRows(uint64_t** pmatrix, int height, int width)
{
	mat input(pmatrix,height, width);
	return input.maxSubmatrix();
}

uint64_t** mat::toNMatrix()
{
	freeMat();
	nMatHeight = getHeight();
	nMat = new uint64_t* [getHeight()];
	for (int i = 0; i < getHeight(); i++)
	{
		nMat[i]=new uint64_t[getWidth()];
	}
	for (int i=0;i<getHeight();i++)
	{
		for (int j = 0; j < getWidth(); j++)
		{
			nMat[i][j] = at(i,j);
		}
	}
	return nMat;
}

void mat::freeMat()
{
	if (nMat != nullptr)
	{
		for (int i = 0; i < nMatHeight; i++)
		{
			delete[] nMat[i];
		}
		delete[] nMat;
	}
}

/*std::pair<std::vector<int>,mat> mat::inverseSubmatrix()
{
	auto uTriangle = upper_triangle_transform();
	std::vector<int> indicies;
	mat triangle;
	mat pInverse;
	std::tie(triangle,pInverse,indicies)= uTriangle;
	std::vector<int> firstNumbers(triangle.getHeight());
	for (int i=0;i<firstNumbers.size();i++)
	{
		firstNumbers[i]=i;
	}
	std::vector<int> remainingNumbers;
	for (int i=triangle.getHeight();i<triangle.getWidth();i++)
	{
		remainingNumbers.push_back(i);
	}
	pInverse.extractMatrix(remainingNumbers,indicies);
	indicies.resize(triangle.getHeight());
	mat inverse = backSubstituation(make_tuple(triangle.extractColumns(firstNumbers),pInverse ,firstNumbers));
	return make_pair(indicies,inverse);
}*/

mat mat::rearrangeMatrix(const std::vector<int> & arrangement)
{
  mat result(getHeight(),getWidth());
  for(int i=0;i<arrangement.size();i++)
  {
    result.col(arrangement[i])=col(i);
  }
  return result;
}
