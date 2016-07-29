#include "matrix.h"
#include "galois.h"

 /**
	* @brief creates a matrix with given height and width from a given array
	*/
mat::mat(uint64_t** pmatrix, int height, int width) 
{
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			addValue(pmatrix[i][j],j);
		}
	}
}

 /**
	* @brief creates a matrix with given height and width
	*/
mat::mat(int height, int width)
{
	for (int i = 0; i<width; i++)
	{
		std::vector<uint64_t>column(height,0);
		matrix.push_back(column);
	}
	updateDimension();
}

 /**
	* @brief creates a matrix
	*/
mat::mat()
{
	updateDimension();
}

 /**
	* @brief creates a matrix
	*/
mat::mat(std::initializer_list<std::initializer_list<uint64_t>> lst)
	: n_rows(lst.size()), n_cols(n_rows ? lst.begin()->size() : 0)
{
	int i = 0;
	for (const auto& row : lst) {
		for (const auto& value : row) {
			addValue(value,i);
			i++;
		}
		i=0;
	}
}

 /**
	* @brief append a value tho the matrix in the given column
	*/

void mat::addValue(uint64_t value, int columnNumber)
{
	if (matrix.size()<= (unsigned) columnNumber)
	{
		matrix.push_back(std::vector<uint64_t>(1, value));
	}
	else
	{
		matrix[columnNumber].push_back(value);
	}
}
 /**
	* @brief zeros the entry of the matrix
	*/
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
 /**
	* @brief returns the entry of the matrix at position a,b
	*/
uint64_t & mat::operator()(const int a,const int b)
{
	return matrix[b][a];
}
 /**
	* @brief returns the entry of the matrix at position a,b
	*/
uint64_t mat::operator()(const int a, const int b) const
{
	return matrix[b][a];
}
 /**
	* @brief returns the entry of the matrix at position a,b
	*/
uint64_t mat::at(const int a, const int b) const
{
	return matrix[b][a];
}
 /**
	* @brief returns the entry of the matrix at position a,b
	*/
uint64_t & mat::at(const int a, const int b)
{
	return matrix[b][a];
}

 /**
	* @brief multiplicates two matrices
	*/
mat mat::operator*(const mat & rhs)
{
	if (getWidth() != rhs.getHeight())
	{
		throw;
	}
	mat result(getHeight(), rhs.getWidth());
	for (std::size_t i = 0; i < result.getWidth(); i++)
	{
		for (std::size_t j=0;j<result.getHeight();j++)
		{
			for (std::size_t k=0;k<getWidth();k++)
			{
				result(j,i) = Galois::getInstance().add(result(j,i),Galois::getInstance().multiply(matrix[k][j],rhs(k,i)));
			}
		}
	}
	return result;
}
 /**
	* @brief adds two matrices
	*/
mat mat::operator+(const mat & rhs)
{
	if (getWidth() != rhs.getWidth() || getHeight()!= rhs.getHeight())
	{
		throw;
	}
	mat result(getHeight(), getWidth());
	for (std::size_t i = 0; i < result.getWidth(); i++)
	{
		for (std::size_t j = 0; j<result.getHeight(); j++)
		{
			result(j,i) = Galois::getInstance().add(matrix[i][j], rhs(j,i));
		}
	}
	return result;
}

 /**
	* @brief subtract two matrices
	*/
mat mat::operator-(const mat & rhs)
{
	return (*this)+rhs;
}
 /**
	* @brief concat the left matrix with the right matrix
	*/
mat join_rows(const mat & left, const mat & right)
{
	if (left.getHeight() != 0 && right.getHeight() != 0 && left.getHeight() != right.getHeight())
	{
		throw;
	}
	mat result=left;
	for (std::size_t i=0;i<right.getWidth();i++)
	{
		result.addColumn(right.col(i));
	}
	return result;
}
/**
	* @brief concat the left matrix a column
	*/
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

 /**
	* @brief updates the saved dimensions of the matrix
	*/
void mat::updateDimension()
{
	n_rows = getHeight();
	n_cols = getWidth();
}
 /**
	* @brief returns a specific column of the matrix
	*/
std::vector<uint64_t> & mat::col(const int colnumber)
{
	return matrix[colnumber];
}
 /**
	* @brief returns a specific column of the matrix
	*/
std::vector<uint64_t> mat::col(const int colnumber) const
{
	return matrix[colnumber];
}
 /**
	* @brief deletes a specific row of the matrix
	*/
void mat::shed_row(const int rowNumber)
{
	for (auto & row : matrix)
	{
		row.erase(row.begin() + rowNumber);
		row.shrink_to_fit();
	}
	updateDimension();
}

 /**
	* @brief adds a specific column
	*/
void mat::addColumn(const std::vector<uint64_t> & column)
{
	matrix.push_back(column);
	updateDimension();
}
 /**
	* @brief prints the matrix
	*/
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
		for (std::size_t j = 0; j<getHeight(); j++)
		{
			output.str(std::string());
			for (std::size_t i = 0; i < getWidth(); i++)
			{
				output<<matrix[i][j] <<" ";
			}
			std::cout<<output.str() << std::endl;
		}
	}
}
 /**
	* @brief calculates the rank of the matrix
	*/
int matRank(mat matrix)
{
	auto transform = matrix.upper_triangle_transform();
	int rank = 0;
	for (std::size_t i = 0; i < std::min(std::get<0>(transform).getWidth(),std::get<0>(transform).getHeight()); i++)
	{
		if (std::get<0>(transform)(i,i)!=0)
		{
			rank++;
		}
	}
	return rank;
}

 /**
	* @brief checks if the matrix is a square matrix
	*/
bool mat::is_square() const
{
	return getWidth() == getHeight();
}

 /**
	* @brief transposes the matrix
	*/
mat mat::t()
{
	mat result(getWidth(), getHeight());
	for (std::size_t i = 0;i<getWidth();i++)
	{
		for (std::size_t j=0;j<getHeight();j++)
		{
			result(i,j) = matrix[i][j];
		}
	}
	return result;
}

 /**
	* @brief invertes the matrix
	*/
mat mat::i()
{
	if (!is_square())
	{
		throw;
	}
	auto uTriangle = upper_triangle_transform();
	return backSubstituation(uTriangle);
}
 /**
	* @brief finds the first non-zero value of a matrix in some row starting right to a specific position
	* @param [in] input the matrix we are interested in
    * @param [in] row the row we want to find the nonzero entry
	* @param [in] the starting colomn number, we will not look left from this column
	* @returns column number if successful, -1 else
	*/
int findNonZero(mat & input,int row, int startCol)
{
	for (std::size_t i = startCol;i<input.getWidth();i++)
	{
		if (input(row, i) != 0)
		{
			return i;
		}
	}
	return -1;
}
 /**
	* @brief finds a maximum submatrix that has full rank
	* @returns indices of rows and columns the submatrix is consisting of
	*/
std::vector<int> mat::maxSubmatrix()
{
	if (getWidth() != getHeight())
	{
		throw std::runtime_error("Wrong Dimensions");;
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

 /**
	* @brief use column operations to zero destColNumber using the pivot element at rowNumber,colNumber. Operate also on inverse
	*/
void mat::columnTransform(mat & matrix, mat & inverse, int colNumber, int destColNumber, int rowNumber)
{
	if (matrix(rowNumber, colNumber) != 1)
	{
		uint64_t inv=Galois::getInstance().inverse(matrix(rowNumber, colNumber));
		columnOperation(matrix,colNumber, inv);
		columnOperation(inverse,colNumber, inv);
	}
	uint64_t fac= matrix(rowNumber,destColNumber);
	columnOperation(matrix, colNumber, destColNumber, fac);
	columnOperation(inverse,colNumber, destColNumber, fac);
}
 /**
	* @brief multiply a column of matrix with a factor
	*/
void mat::columnOperation(mat& matrix, int colNumber, uint64_t factor)
{
	for (std::size_t i = 0; i < matrix.getHeight(); i++)
	{
		matrix(i,colNumber) = Galois::getInstance().multiply(matrix(i,colNumber),factor);
	}
}
 /**
	* @brief multiply the first column with a factor and add to the second column
	*/
void mat::columnOperation(mat & matrix, int firstColNumber, int secondColNumber, uint64_t factor)
{
	for (std::size_t i = 0; i < matrix.getHeight(); i++)
	{
		matrix(i,secondColNumber) = Galois::getInstance().add(Galois::getInstance().multiply(matrix(i,firstColNumber), factor), matrix(i,secondColNumber));
	}
}
 /**
	* @brief use row operations to zero destRowNumber using the pivot element at rowNumber,colNumber. Operate also on inverse
	*/
void mat::rowTransform(mat & matrix, mat & inverse, int rowNumber, int destRowNumber, int colNumber)
{
	if (matrix(rowNumber, colNumber) != 1)
	{
		uint64_t inv = Galois::getInstance().inverse(matrix(rowNumber, colNumber));
		rowOperation(matrix, rowNumber, inv);
		rowOperation(inverse, rowNumber, inv);
	}
	uint64_t fac = matrix(destRowNumber, colNumber);
	rowOperation(matrix, rowNumber, destRowNumber, fac);
	rowOperation(inverse, rowNumber, destRowNumber, fac);
}
 /**
	* @brief use row operations to zero destRow using the pivot element at rowNumber,colNumber
	*/
void mat::rowTransform(mat & matrix, int rowNumber, int destRowNumber, int colNumber)
{
	if (matrix(rowNumber, colNumber) != 1)
	{
		uint64_t inv = Galois::getInstance().inverse(matrix(rowNumber, colNumber));
		rowOperation(matrix, rowNumber, inv);
	}
	uint64_t fac = matrix(destRowNumber, colNumber);
	rowOperation(matrix, rowNumber, destRowNumber, fac);
}
 /**
	* @brief multiply a row of matrix with a factor
	*/
void mat::rowOperation(mat& matrix, int rowNumber, uint64_t factor)
{
	for (size_t i = 0; i < matrix.getWidth(); i++)
	{
		if (factor == 0)
		{
			matrix(rowNumber, i) = 0;
		}
		else if (factor!=1)
		{
			matrix(rowNumber, i) = Galois::getInstance().multiply(matrix(rowNumber, i), factor);
		}
	}
}
 /**
	* @brief multiply the frist row with a factor and add to the second row
	*/
void mat::rowOperation(mat & matrix, int firstRowNumber, int secondRowNumber, uint64_t factor)
{
	for (std::size_t i = 0; i < matrix.getWidth(); i++)
	{
		if (factor == 1)
		{
			matrix(secondRowNumber, i) = Galois::getInstance().add(matrix(firstRowNumber, i), matrix(secondRowNumber, i));
		}
		else if (factor != 0)
		{
			matrix(secondRowNumber, i) = Galois::getInstance().add(Galois::getInstance().multiply(matrix(firstRowNumber, i), factor), matrix(secondRowNumber, i));
		}
	}
}

 /**
	* @brief calculates the determinant of the matrix
	*/

uint64_t mat::det()
{
	auto transform = upper_triangle_transform();
	uint64_t result = std::get<0>(transform)(0,0);
	for (std::size_t i = 1; i < std::get<0>(transform).getWidth(); i++)
	{
		result = Galois::getInstance().multiply(result, std::get<0>(transform)(i,i));
	}
	return result;
}

 /**
	* @brief swaps the given two rows
	*/

void mat::swapColumns(const int a,const int b)
{
	auto temp = matrix[a];
	matrix[a] = matrix[b];
	matrix[b]=temp;
}

 /**
	* @brief backSubstituation step of the Gauss algorithm
	* @param [in] upperTriangle the matrix, inverted matrix so far and permutation of the columns
	* @returns the inverse
	*/
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

 /**
	* @brief transform matrix to upper triangle form
	* @returns matrix, inverted matrix so far, permutation of the columns
	*/
std::tuple<mat,mat,std::vector<int>> mat::upper_triangle_transform()
{
	mat matrix = (*this);
	mat pInverse = eye<mat>(getWidth(), getWidth());
	std::vector<int> result;
	for (std::size_t i = 0; i < getWidth(); i++)
	{
		result.push_back(i);
	}
	for (std::size_t n = 0; n < std::min(matrix.getHeight(),matrix.getWidth()); n++)
	{
		//find a nonzero entry
		int nonzero = findNonZero(matrix, n, n);
		if (nonzero!=-1)
		{
			//swap the columns such that 
			if ((unsigned) nonzero != n)
			{
				matrix.swapColumns(nonzero, n);
				pInverse.swapColumns(nonzero, n);
				std::swap(result[nonzero],result[n]);
			}
			//use column transformation to zero entries right to the diagonal
			for (std::size_t i = n + 1; i < matrix.getWidth(); i++)
			{
				columnTransform(matrix, pInverse, n, i, n);
			}
		}
	}
	return std::make_tuple(matrix,pInverse,result);
}

 /**
	* @brief transform the matrix to the from (IB) where I is the identity matrix
	* @returns the matrix and the permutation of the columns
	*/
std::pair<mat,std::vector<int>> mat::toStandarForm()
{
	mat matrix = (*this);
	std::vector<int> result;
	for (std::size_t i = 0; i < getWidth(); i++)
	{
		result.push_back(i);
	}
	for (std::size_t n = 0; n < std::min(matrix.getHeight(), matrix.getWidth()); n++)
	{
		//find nonzero entry
		int nonzero = findNonZero(matrix, n, n);
		if (nonzero != -1)
		{
			//swap columns if necessary
			if ((unsigned) nonzero != n)
			{
				matrix.swapColumns(nonzero, n);
				std::swap(result[nonzero], result[n]);
			}
			//use row operation to zero entries right the diagonal
			for (std::size_t i = n + 1; i < matrix.getHeight(); i++)
			{
				rowTransform(matrix, n, i, n);
			}
		}
	}
	//use row transformation to zero entries not on the diagonal	
	for (int i = matrix.getHeight() - 1; i >= 0; i--)
	{
		for (int j = 0; j<i; j++)
		{
			rowTransform(matrix, i, j, i);
		}
	}
	return std::make_pair(matrix, result);
}

 /**
	* @brief matracts the given columns of the matrix
	*/
mat mat::extractColumns(const std::vector<int> & index)
{
	mat result;
	for (int i : index)
	{
		result = join_rows(result, col(i));
	}
	return result;
}

 /**
	* @brief compares two matrices
	*/
bool mat::operator==(const mat & rhs) 
{
	if (getHeight() != rhs.getHeight() || getWidth() != rhs.getWidth())
	{
		return false;
	}
	for (std::size_t i=0;i<getHeight();i++)
	{
		for (std::size_t j=0;j<getWidth();j++)
		{
			if (matrix[j][i] != rhs(i,j))
			{
				return false;
			}
		}
	}
	return true;
}

 /**
	* @brief compares to matrices
	*/
bool mat::operator!=(const mat & rhs)
{
	return !((*this) == rhs);
}

 /**
	* @brief generates a matrix from array and return indices a maximal full rank submatrix is consisting of
	*/
std::vector<int> findRedundantRows(uint64_t** pmatrix, int height, int width)
{
	mat input(pmatrix,height, width);
	return input.maxSubmatrix();
}

 /**
	* @brief converts the matrix to array
	*/
uint64_t** mat::toNMatrix()
{
	uint64_t** nMat = new uint64_t* [getHeight()];
	for (std::size_t i = 0; i < getHeight(); i++)
	{
		nMat[i]=new uint64_t[getWidth()];
	}
	for (std::size_t i=0;i<getHeight();i++)
	{
		for (std::size_t j = 0; j < getWidth(); j++)
		{
			nMat[i][j] = at(i,j);
		}
	}
	return nMat;
}

 /**
	* @brief rearranges the columns of the matrix by the given permutation
	*/
mat mat::rearrangeMatrix(const std::vector<int> & arrangement)
{
  mat result(getHeight(),getWidth());
  for(std::size_t i=0;i<arrangement.size();i++)
  {
    result.col(arrangement[i])=col(i);
  }
  return result;
}
