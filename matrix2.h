#ifndef MATRIX_H_
#define MATRIX_H_


#include <iostream>

#include <cstdlib> 
#include <ctime>

/**
	* @brief This function multiplies a scalar to a matrix.
	*
	* Manipulates the delivered matrix by multiplying all values with a scalar.
	*
	* @param [in] gal The Galois field structure.
	* @param [in] mat matrix of size row cross col
    * @param [in] scalar the scalar the matrix is multiplied with
    * @param [in] row the number of rows of the delivered matrix
    * @param [in] col the number of collumns of the delivered matrix
	*/
void scalar_matrix(Galois gal, uint64_t** mat, uint64_t scalar, int row, int col);

/**
	* @brief This function add a matrix to a matrix.
	*
	* Manipulates the first delivered matrix by adding the second matrix.
	*
	* @param [in] gal The Galois field structure.
	* @param [in] mat1 matrix of size row x col, which will be manipulated by adding the second matrix
    * @param [in] mat2 matrix of size row x col
    * @param [in] row the number of rows of the delivered matrix
    * @param [in] col the number of collumns of the delivered matrix
	*/
void add_matrix(Galois gal, uint64_t** mat1, uint64_t** mat2, int row, int col);


/**
	* @brief Prints the delivered matrix to the console.
	*
	* 
	*
	* @param [in] gal The Galois field structure.
    * @param [in] size_row the number of rows of the delivered matrix
    * @param [in] size_col the number of collumns of the delivered matrix
    * @param [in] matrix matrix of size row x col
	*/
void print_matrix(Galois g, int size_row, int size_col, uint64_t** matrix);



/**
	* @brief Prints the delivered vector (uint64_t) to the console.
	*
	* 
	*
	* @param [in] gal The Galois field structure.
    * @param [in] size the size of the vector 
    * @param [in] vector vector of size size
	*/
void print_vector(Galois g, int size, uint64_t* vector);

/**
	* @brief Prints the delivered vector (int) to the console.
	*
	* 
	*
    * @param [in] size the size of the vector 
    * @param [in] vector vector of size size
	*/
void print_vector_normal(int size, int* vector);


/**
	* @brief Swaps two rows of the matrix
	*
	* Manipulates the delivered matrix by swaping row "line1" with the row "line2".
	* If the indixes are not in matrix, it will return false, otherwise (successfully swapped) true
	*
	* @param [in] mat The matrix of size row x col
	* @param [in] row the number of rows of the delivered matrix
    * @param [in] col the number of collumns of the delivered matrix
    * @param [in] line1 the index of the first row to swap
    * @param [in] line2 the index of the second row to swap
	* @returns true if sussessfull, false otherwise
	*/
bool swapLine(uint64_t** mat, int row, int col, int line1, int line2);

 /**
	* @brief Invert a matrix.
	*
	* Computes the inverted matrix of a size x size matrix by using the Gauß-Jordan-Algorithm.
	*
	* @param [in] gal The Galois field structure.
	* @param [in] mat matrix of size size x size
    * @param [in] size the size of the squared matrix 
	* @returns the inverted matrix
	*/
uint64_t** invertMatrix(Galois gal, uint64_t** mat, int size);


 /**
	* @brief Multiplies two matrices 
	*
	* Computes the product of two matrices by creating a new one.
	*
	* @param [in] gal The Galois field structure.
	* @param [in] A matrix of size row_A x col_A
	* @param [in] row_A number of rows of the first matrix
	* @param [in] col_A number of collumns of the first matrix
	* @param [in] B matrix of size row_A x col_A
	* @param [in] row_B number of rows of the second matrix
	* @param [in] col_B number of collumns of the second matrix

	* @returns the row_A x col_B sized matrix. which is the product of A times B
	*/
uint64_t** multiplication_matrix(Galois gal, uint64_t** A, int row_A, int col_A, uint64_t** B, int row_B, int col_B);



  /**
	* @brief An auxiliary function for the small-rank-update formula
	*
	* Computes I + V_transposed * M * U with V is size 2 x n and U is size n x 2 and M is size n x n.
	* This is needed in the small rank update formula by Sherman,Morrison and Woodbury
	*
	* @param [in] gal The Galois field structure.
	* @param [in] V matrix of size 2 x size
	* @param [in] M matrix of size size x size
	* @param [in] U matrix of size size x 2
	* @param [in] size the needed parameter for the matrix structure


	* @returns I + V^T M U
	*/
uint64_t** SMW_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** U, int size);


   /**
	* @brief An auxiliary function for the small-rank-update formula
	*
	* Computes M U (I + V M U) V M  with V is size 2 x n and U is size n x 2 and M is size n x n.
	* This is needed in the small rank update formula by Sherman,Morrison and Woodbury
	*
	* @param [in] gal The Galois field structure.
	* @param [in] V matrix of size 2 x size
	* @param [in] M matrix of size size x size
	* @param [in] SMW matrix of size size x size (the result of SMW_matrix method)
	* @param [in] U matrix of size size x 2
	* @param [in] size the needed parameter for the matrix structure


	* @returns  M U (I + V M U) V M  
	*/
uint64_t** SMW_inverse_update_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** SMW, uint64_t** U, int size );

   /**
	* @brief Computes the wedge product of two vectors
	*
	* Computes the wedge product of b and c, which are vector of size size
	*
	* @param [in] gal The Galois field structure.
	* @param [in] b vector of size size
	* @param [in] c vector of size size
	* @param [in] size the needed parameter for the vector structure


	* @returns  size x size matrix, which is the wedge product of the vectors b and c
	*/
uint64_t** wedge_product(Galois gal, uint64_t *b, uint64_t *c, int size);

   /**
	* @brief Get a copy of a matrix
	*
	* Allocates and assign a new matrix with the same values as the delivered one
	*
	*
	* @param [in] mat The matrix of size row x col
	* @param [in] row the number of rows of the delivered matrix
    * @param [in] col the number of collumns of the delivered matrix


	* @returns  copy of the delivered matrix
	*/
uint64_t** copy_matrix(uint64_t** mat, int row, int col);

   /**
	* @brief Frees the allocated memory
	*
	* Frees the allocated memory of the matrix
	*
	* @param [in] mat The matrix with numer of rows is "row"
	* @param [in] row the number of rows of the delivered matrix
  
	*/
void my_free(uint64_t** mat, int row);


   /**
	* @brief compares two matrices component for component
	*
	* Checks, if the two given matrices are equal. (Used for debug)
	*
	* @param [in] gal The Galois field structure.
	* @param [in] A The first matrix of size "size" 
	* @param [in] B The secound matrix of size "size"
	* @param [in] size the size of the delivered matrices
  	*
  	* @returns true if same, false if unequal
	*/
bool compare(Galois gal, uint64_t** A, uint64_t** B, int size);
#endif /* LIN_PARITY_H_ */