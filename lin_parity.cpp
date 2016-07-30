
#include "galois.h"
#include "gauss.h"
#include "lin_parity.h"
#include "matrix.h"
#include "matrix2.h"

/**
	* @brief Creates a random value
	* Creates a non-zero integral random value between 1 and max.
	* @param [in] max THe upper bound for the integral random value.
	* @return random value
	*/
int get_random_value(int max){
	int random_variable = std::rand();
    return random_variable % max + 1 ;
}

/* Creates the compact Matrix Y and stores the random values x_i in vector random_values
 * so random_values has to be a "column of M-half-dim" vec (a pointer of it)
 * random values are integral between 1 and max_random_value
 */ 

 /**
	* @brief Create the compact matrix representation for the linear matriod parity problem
	*
	* Creates matrix Y of algorithm 4.1, the compact representation of the linear matroid parity problem.
	* The used random values will be stored in the vector "random_values" (for later use).
	* the random value
	*
	* @param [in] mat The matrix of size row x col
	* @param [in] row the number of rows of the delivered matrix
    * @param [in] col the number of collumns of the delivered matrix
    * @param [in] line1 the index of the first row to swap
    * @param [in] line2 the index of the second row to swap
	* @returns true if sussessfull, false otherwise
	*/
uint64_t** create_Y(Galois & gal, uint64_t **M, int row, int col, uint64_t *random_values){
	
	//create a new row x row matrix (will be Y, the result)
	uint64_t** result_matrix= new uint64_t*[row];
	for (int i=0; i<row; i++){
    	result_matrix[i]= new uint64_t[row];
	}
	
	//initialise the result_matrix with zeros
	for (int i=0; i<row; i++){
    	for(int j=0; j<row; j++){
        	result_matrix[i][j]= (uint64_t) 0;
    	}
	}
	
	for(int i = 0 ; i < col ; i = i+2){
		///compute Y = Y + random_value (b ^ c) where b is the i-th col and c the i+1-th col
		while (random_values[i/2] == 0 ){ 
			
			random_values[i/2] = gal.uniform_random_element(); // if no specific random_value vector is handed, a new random value is assigned 	(!=0)
			
		}
		

		for(int row_index = 0 ; row_index < row; row_index++){		//row is dimension of Y (Y is row x row matrix, want to update Y)
			for (int col_index = 0; col_index < row ; col_index++){
				//want to add x_i (M[][i] ^ M[][i+1]) to Y
				result_matrix[row_index][col_index] =gal.add(result_matrix[row_index][col_index], gal.multiply(random_values[i/2], gal.add(gal.multiply(M[row_index][i], M[col_index][i+1]), gal.multiply(M[row_index][i+1], M[col_index][i] ) )  ));
			}
		}
	
	}

	return result_matrix;
}

/**
	* @brief Creates the matrix U for the algo 4.1.
	*
	* Creates the matrix U for the algo 4.1, so it takes the i-th and the (i+1)-th column of M, scaled by the random value
	* U = x_i * ( i-th col  |  i+1-th col )
	*
	* @param [in] gal The galois field structure
	* @param [in] M The main matrix of lin parity 
    * @param [in] size the number of rows of M
    * @param [in] random the random value, the result size x 2 matrix will be scaled with this 
    * @param [in] i the number of the desired i-th collumn
	* @returns The Matrix U
	*/
uint64_t** get_U(Galois & gal, uint64_t** M, int size, uint64_t random, int i){
	//construce a size x 2 matrix
	uint64_t** U = new uint64_t*[size];
	for (int k=0; k<size; k++){
    	U[k]= new uint64_t[2];
	}

	for(int k = 0; k < size; k++){
		U[k][0] = M[k][i];
		U[k][1] = M[k][i+1];
	}
	
	scalar_matrix(gal,  U, random, size, 2);
	return U;
}

/**
	* @brief Creates the matrix V for the algo 4.1.
	*
	* Creates the matrix V for the algo 4.1, so it takes the i+1-th and the i-th column of M and uses the transposed
	* V = ( i+1-th col   i-th col) ^T
	*
	* @param [in] gal The galois field structure
	* @param [in] M The main matrix of lin parity 
    * @param [in] size the number of rows of M 
    * @param [in] i the number of the desired i-th collumn
	* @returns The Matrix V
	*/
 uint64_t** get_V(Galois & gal, uint64_t** M, int size, int i){

 	//construce a 2 x size matrix
	uint64_t** V = new uint64_t*[2];
	for (int k=0; k<2; k++){
    	V[k]= new uint64_t[size];
	}
	
	for (int k = 0; k <size ; k++){
		V[0][k] = M[k][i+1];
		V[1][k] = M[k][i];
	}
	
	return V;
 }



/*
 * This implemetation uses the small rank update formula and also return a maximum number of pairs, if there is no
 * parity basis
 * in length it returns the number of lin. independent vectors (so there are lenght/2 pairs)
 */
 /**
	* @brief Computes the max number of linear independent pairs of the delivered matrix M
	*
	* Solves the linear matroid parity problem for the given instance M, so it detects the maximum number of pairs, which are
	* lineary independent. It uses the small rank update formula by sherman-morrison-woodbury to not compute in every itaration 
	* the determinant from scratch.
	*
	* @param [in] gal The galois field structure
	* @param [in] M The main matrix of lin parity 
    * @param [in] row the number of rows of M
    * @param [in] col the number of rows of M 
    * @param [in] length here we will store the length of the output (number of lin indepentent vectors )
	* @returns an array with all the numbers of cols of the lin independent vector-pairs
	*/
int* simple_parity_fast(Galois & gal, uint64_t** M, int row, int col, int* length){
	int success = 0;
	int* parity_basis;

while (success == 0){
	success = 1;

	int check = 0;
	int row_difference = 0;
	Gauss gauss;
	uint64_t* random_values = new uint64_t[col/2]; 		//stores the randomvalues x_i, created in create_Y (needed later on) 
	for (int i = 0; i < col/2; i++){
		random_values[i] = 0;
	}
	
	

	/*cout<<"------------------Algorithm start-----------"<< endl;
	cout << "Startmatrix M :"<<endl;
	if (row < 13) {debug print_matrix(gal,row, col, M);}
	else {cout << "too big, will not write it down" <<endl;}
	*/
	//computing Y
	uint64_t** Y = create_Y( gal, M, row, col, random_values);
	uint64_t** Y_copy_det = copy_matrix(Y, row, row); //copy matrix for computing the determinant (matrix will get transformed by the det function)

	uint64_t det = gauss.determinant(gal, row, Y_copy_det);	
	
	my_free(Y_copy_det, row);
	
	if (det == 0) {
		check = 1;
		// There is no parity basis, so we are searching for redundent rows to delete
		
		std::vector<int> del_row  = findRedundantRows(Y, row, row);
		
		row -= del_row.size();   				//updating our new number of rows
		row_difference = del_row.size();		//just needed for correctly freeing space later	
		sort(del_row.begin(), del_row.end());
		
		if (row == 0){
			//everything was redunden, new matrix is empty, return this
			delete[] random_values;
			my_free(Y, row + row_difference);
			my_free(M, row + row_difference);
			int* parity_basis = new int[row];		    //in here we will store the indices of the columns that built the parity basis
			*length = 0;
			return parity_basis;
		}


		//delete the rows in M
		uint64_t** M_prime = new uint64_t*[row];
			for (int k=0; k<row; k++){
    			M_prime[k]= new uint64_t[col];
		}

		unsigned int counter = 0;
		for (unsigned int i = 0; i < row + del_row.size() ; i++ ){
			if ((counter < del_row.size()) &&   (i == (unsigned) del_row[counter])) {
				counter++;
			}
			else{
				for(int j = 0; j < col; j++){			
					M_prime[i-counter][j] = M[i][j];	
				}
			}
		}

		my_free(M, row + row_difference);

		//now we have our new (redundent-free, full rank) matrix M with det(Y) != 0 
		M = M_prime;

	}


	//have to compute Y once again for the new M
	if (check == 1){
		my_free(Y, row + row_difference);
		for (int i = 0; i < col/2; i++){
			random_values[i] = 0;
		}
		Y = create_Y( gal, M, row, col, random_values);
	}

	uint64_t** Y_inverse = invertMatrix(gal, Y,  row);

	if (Y_inverse==0)
	{
		success=0;
		my_free(Y, row);
		continue;
	}

	parity_basis = new int[row];		    //in here we will store the indices of the columns that built the parity basis
	for (int i = 0; i < row; i++){
		parity_basis[i] = 0;
	}
	
	int pos = 0;			   						//for positioning vector_parity 
	for(int i = 0; i < col; i = i+2){
		//using the small rank update formula (instead of computing det of Y_prime we compute it of a 2x2 matrix!)
		uint64_t** V = get_V(gal, M, row, i);
		uint64_t** U = get_U(gal, M, row,random_values[i/2], i);	

 		uint64_t** SMW_det = SMW_matrix(gal,  V, Y_inverse, U, row);
		uint64_t** SMW = copy_matrix(SMW_det, 2,2);
		
		uint64_t det_SMW = gauss.determinant(gal, 2, SMW_det);
		
		my_free(SMW_det,2);

		if (det_SMW != 0){ //if det(Y_prime) != 0

			//Setzt Y auf Y' ( also Y + UV)
			uint64_t** hilf = multiplication_matrix(gal, U, row, 2, V, 2, row);
			add_matrix(gal, Y, hilf, row, row);

			//update Y_inverse
			uint64_t** update_matrix = SMW_inverse_update_matrix( gal,  V,  Y_inverse, SMW,  U, row );
			add_matrix( gal, Y_inverse, update_matrix,  row, row);
			my_free(U, row);
			my_free(V, 2);
			my_free(SMW, 2);
			my_free(hilf, row);
			my_free(update_matrix, row);
		
		}
		else{
			my_free(U, row);
			my_free(V, 2);
			my_free(SMW, 2);
			if (pos>=row)
			{	
				success = 0;
				break;
			}
			parity_basis[pos] = i;
			parity_basis[pos+1] = i+1;
			pos = pos + 2;
		}
	}
	delete [] random_values;
	*length = row;

	my_free(Y, row);
	my_free(Y_inverse, row);


	if (success == 1) my_free(M, row);

	
 }

 return parity_basis;
}







/*int main(){
	
	Galois gal;
  	Gauss gauss;
  	gal.set_w(16);
  
  	gal.set_mode_logtb();
  
  	gal.seed();
	
	int row = 41;
	int col = 60;
	
 
 	//construct a row x col matrix
	uint64_t** matrix= new uint64_t*[row];
 
	for (int i=0; i<row; i++){
    	matrix[i]= new uint64_t[col];
	}
 	
	//fill it with random numbers
	for (int i=0; i<row; i++){
    	for(int j=0; j<col; j++){
        	matrix[i][j]= gal.uniform_random_element();
        	if (j == (col-1)){
        		 matrix[i][j-3] = matrix[i][j];
        		 matrix[i][j-2] = matrix[i][j];
        		 matrix[i][j-1] = matrix[i][j];
        		
        	}
    	}
	}
	cout<< "Startmatrix = ";
	print_matrix(gal, row, col, matrix);
	
	int length;

	int* result = simple_parity_fast(gal, matrix, row, col, &length);

	cout << "ENDERGEBNIS " << endl;
	print_vector_normal(length, result);

  	my_free(matrix, row);	
  	delete [] result;
	return 0;

}*/