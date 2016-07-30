#include "galois.h"
#include "gauss.h"
#include "lin_parity.h"
#include "matrix2.h"


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
void scalar_matrix(Galois gal, uint64_t** mat, uint64_t scalar, int row, int col){
	for(int i = 0; i < row ; i++){
		for(int j = 0; j < col; j++){
			mat[i][j] = gal.multiply(mat[i][j], scalar);
		}
	}

}

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
void add_matrix(Galois gal, uint64_t** mat1, uint64_t** mat2, int row, int col){
	for(int i = 0; i < row ; i++){
		for(int j = 0; j < col; j++){
			mat1[i][j] = gal.add(mat1[i][j], mat2[i][j]);
		}
	}
	
}

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
void print_matrix(Galois g, int size_row, int size_col, uint64_t** matrix)
{
  cout<<"Matrix is:" <<endl;

  for (int i = 0; i < size_row; i++){
  	
  	for (int j = 0; j < size_col; j++){
  	
  		cout << g.to_string(matrix[i][j]) << " " ;
  	}
  	cout << endl;
  	
  }
  
  
}



/**
  * @brief Prints the delivered vector (uint64_t) to the console.
  *
  * 
  *
  * @param [in] gal The Galois field structure.
    * @param [in] size the size of the vector 
    * @param [in] vector vector of size size
  */
void print_vector(Galois g, int size, uint64_t* vector)
{
  cout<<"Vector is:" <<endl;
  for (int i = 0; i < size; i++){
  	
  		cout << g.to_string(vector[i]) << endl;
 
  }
  
}

/**
  * @brief Prints the delivered vector (int) to the console.
  *
  * 
  *
    * @param [in] size the size of the vector 
    * @param [in] vector vector of size size
  */
void print_vector_normal(int size, int* vector)
{
  cout<<"Vector is:" <<endl;
  for (int i = 0; i < size; i++){
  	
  		cout <<(vector[i]) << " ";
   
  }
  cout <<endl;
  
}

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
bool swapLine(uint64_t** mat, int row, int col, int line1, int line2)
{
    if(line1 >= row || line2 >= row)
        return false;
 
    for(int i = 0; i < col; ++i)
    {
        uint64_t t = mat[line1][i];
        mat[line1][i] = mat[line2][i];
        mat[line2][i] = t;
    }
     
    return true;
}



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
uint64_t** invertMatrix(Galois gal, uint64_t** mat, int size)
{

  //get space for the inverted matrix
	uint64_t** inv= new uint64_t*[size];
	for (int i=0; i<size; i++){
    	inv[i]= new uint64_t[size];
	}
   
   //Get a size x 2*size matrix for the Gauß-Jordan_Algorithm
    uint64_t** A= new uint64_t*[size];
	for (int i=0; i<size; i++){
    	A[i]= new uint64_t[size*2];
	}

    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j)
            A[i][j] = mat[i][j];
        for(int j = size; j < 2*size; ++j)
            A[i][j] = (i==j-size) ? (uint64_t) 1 : (uint64_t) 0;
    }
 
    // Gauß-Algorithm
    for(int k = 0; k < size-1; ++k)
    {
        
      //Swapp rows, if pivot element is zero
        if(A[k][k] == 0)
        {
            for(int i = k+1; i < size; ++i)
            {
                if(A[i][k] != 0)
                {
                    
                    swapLine(A, size, 2*size, k, i);
                    break;
                }
                else if(i==size-1)
                    return 0; // There is no element != 0
            }
        }
 
        
        //Eliminate elements under the pivot element
        for(int i = k+1; i < size; ++i)
        {
            uint64_t p =  gal.divide(A[i][k],A[k][k]);     
            for(int j = k; j < 2*size; ++j)
           
            	A[i][j] = gal.add(A[i][j], gal.multiply(A[k][j], p) );
        }
    }
 
    //Compute determinant
   
    /*uint64_t det = A[0][0];
  	for (int i = 1; i < size; i++)
  	{
    det = gal.multiply(det, A[i][i]);
  	}*/
    bool det=true;
    for (int i = 0; i < size; i++)
  	{
      if (A[i][i]==0)
      {
        det=false;
      }
  	}

    if(det == 0)  // Determinant is 0  -> matrix not invertable
	{
       	my_free(inv, size);
       	my_free(A, size);
        return 0;
	}
 
    // Jordan-part of the algorithm
    for(int k = size-1; k > 0; --k)
    {
        for(int i = k-1; i >= 0; --i)
        {
            uint64_t p = gal.divide(A[i][k],A[k][k]);
            for(int j = k; j < 2*size; ++j)
                
            	A[i][j] = gal.add(A[i][j], gal.multiply(A[k][j], p));
        }
    }
 
    // norm entries of left matrix to one and write them into matrix inv
    for(int i = 0; i < size; ++i)
    {
        uint64_t f = A[i][i];
        for(int j = size; j < 2*size; ++j)
            inv[i][j-size] = gal.divide(A[i][j],f);
    }
 	my_free(A, size);
    return inv;
}

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
uint64_t** multiplication_matrix(Galois gal, uint64_t** A, int row_A, int col_A, uint64_t** B, int row_B, int col_B){
	//construct a row_A x col_B matrix
	uint64_t** C = new uint64_t*[row_A];
	for (int i=0; i<row_A; i++){
    	C[i]= new uint64_t[col_B];
	}

	for(int i = 0; i < row_A; i++){
		for (int j = 0; j < col_B; j++){
			//berechne C[i][j]
			C[i][j] = (uint64_t) 0; //initialise with 0
			//now compute the i-th row of A times the j-th column of B
			for(int k = 0; k < col_A; k++){
				C[i][j] = gal.add(C[i][j], gal.multiply(A[i][k] , B[k][j]) ) ; 
			}

		}
	}

	return C;
}


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
uint64_t** SMW_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** U, int size){


	uint64_t** VM_hilf = multiplication_matrix(gal, V, 2, size, M, size, size);    // 2 x size

	uint64_t** result = multiplication_matrix(gal, VM_hilf, 2, size, U, size, 2);  //size x size

	my_free(VM_hilf, 2);

	//now add 1 to the diagonal
	result[0][0] = gal.add( result[0][0], (uint64_t) 1);
	result[1][1] = gal.add( result[1][1], (uint64_t) 1);
	return result;
}

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
uint64_t** SMW_inverse_update_matrix(Galois gal, uint64_t** V, uint64_t** Y_inverse, uint64_t** SMW, uint64_t** U, int size ){

	uint64_t** MU = multiplication_matrix(gal, Y_inverse, size, size, U, size, 2); //size x 2

	uint64_t** SMW_invert = invertMatrix(gal, SMW, 2);	//2x2		

	uint64_t** VM = multiplication_matrix(gal, V, 2, size, Y_inverse, size, size); // 2 x size

	uint64_t** MU_SMW_invert = multiplication_matrix(gal, MU, size, 2, SMW_invert, 2, 2); //size x 2
	uint64_t** result = multiplication_matrix(gal, MU_SMW_invert, size, 2, VM, 2, size); //size x size

	my_free(MU, size);

	my_free(VM, 2);
	my_free(MU_SMW_invert, size);
	my_free(SMW_invert,2);

	return result;
}

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
uint64_t** wedge_product(Galois gal, uint64_t *b, uint64_t *c, int size){
	//construct a size x size matrix
	uint64_t** bct= new uint64_t*[size];
	for (int i=0; i<size; i++){
    	bct[i]= new uint64_t[size];
	}
	//fill it
	for (int i=0; i<size; i++){
    	for(int j=0; j<size; j++){
        	
        	bct[i][j]= gal.multiply(b[i], c[j]);
 	
    	}
	}

	uint64_t** cbt= new uint64_t*[size];
	for (int i=0; i<size; i++){
    	cbt[i]= new uint64_t[size];
	}
	//fill it
	for (int i=0; i<size; i++){
    	for(int j=0; j<size; j++){
        	
        	cbt[i][j]= gal.multiply(c[i], b[j]);
 	
    	}
	}
	add_matrix(gal,bct, cbt, size, size);
	my_free(cbt, size);
	return bct;

}




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
uint64_t** copy_matrix(uint64_t** mat, int row, int col){
	//construct a row x col matrix
	uint64_t** copy= new uint64_t*[row];
	for (int i=0; i<row; i++){
    	copy[i]= new uint64_t[col];
	}
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			copy[i][j] = mat[i][j];
		}
	}
	return copy;
}

  /**
  * @brief Frees the allocated memory
  *
  * Frees the allocated memory of the matrix
  *
  * @param [in] mat The matrix with numer of rows is "row"
  * @param [in] row the number of rows of the delivered matrix
  
  */
void my_free(uint64_t** mat, int row){
	for (int i=0; i<row; i++)
    	delete [] mat[i];
	delete [] mat;
}


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
bool compare(Galois gal, uint64_t** A, uint64_t** B, int size){
  int i , j;
  bool check = true;
  for ( i =0; i < size; i++){
    for( j = 0; j < size; j++){
      if (gal.add(A[i][j], B[i][j])) check = false;
    }
  }

  return check;
}