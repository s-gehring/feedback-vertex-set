
#include "galois.h"
#include "gauss.h"
#include "lin_parity.h"
/* NOT NEEDED
 /* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
/*
 * Not needed
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}
/*
 * Not needed
int determinant_sign(const permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}
 /*
 * Not needed
 
double determinant( matrix<double>& m ) {
    permutation_matrix<std::size_t> pm(m.size1());
    double det = 1.0;
    if(lu_factorize(m,pm) ) {
        det = 0.0;
    } else {
        for(int i = 0; i < m.size1(); i++) 
            det *= m(i,i); // multiply by elements on diagonal
        det = det * determinant_sign( pm );
    }
    return det;
}

*/


int get_random_value(int max){
	int random_variable = std::rand();
    return random_variable % max + 1 ;
}



/* Creates the compact Matrix Y and stores the random values x_i in vector random_values
 * so random_values has to be a "column of M-half-dim" vec (a pointer of it)
 * random values are integral between 1 and max_random_value
 */ 
 
//uint64_t** create_Y(mat M, vec *random_values, int max_random_value){
uint64_t** create_Y(Galois gal, uint64_t **M, int row, int col, uint64_t *random_values){
	
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
		//random_values[i/2] = (uint64_t) 1;
		cout << "Random value ist " << gal.to_string(random_values[i/2]) << endl;	//if something goes wrong, you can see the random value on the console

		for(int row_index = 0 ; row_index < row; row_index++){		//row ist die dimension von Y (Y is row x row matrix, want to update Y)
			for (int col_index = 0; col_index < row ; col_index++){
				//want to add x_i (M[][i] ^ M[][i+1]) to Y
				result_matrix[row_index][col_index] =gal.add(result_matrix[row_index][col_index], gal.multiply(random_values[i/2], gal.add(gal.multiply(M[row_index][i], M[col_index][i+1]), gal.multiply(M[row_index][i+1], M[col_index][i] ) )  ));
			}
		}

		//cout << endl << "Zwischenergebnis für i = " << i << endl;
		//print_matrix(gal, row, row, result_matrix);

		
	}

	return result_matrix;
}


/*
 * basically not needed, this is a computation with 2 intermediate steps (but this means there a two row x row matrices)
 * in the other computation (ref create_Y() ) there is a computation for every cell in one line
 */
uint64_t** create_Y_naive(Galois gal, uint64_t **M, int row, int col, uint64_t *random_values){

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
		//compute Y = Y + random_value (b ^ c) where b is the i-th col and c the i+1-th col
		while (random_values[i/2] == 0 ) {
			random_values[i/2] = gal.uniform_random_element(); // if no specific random_value vector is handed, a new random value is assigned (!=0)
		}
		random_values[i/2] = (uint64_t) 1;
		cout << "Random value ist " << gal.to_string(random_values[i/2]) << endl;	//if something goes wrong, you can see the random value on the console

		uint64_t* b = new uint64_t[row];
		
		for (int j = 0; j < row; j++){
			b[j] = M[j][i];
		}

		uint64_t* c = new uint64_t[row];
		for (int j = 0; j < row; j++){
			c[j] = M[j][i+1];
		}

		cout << "b ist " << endl;
		print_vector(gal, row, b );
		cout << "c ist " << endl;
		print_vector(gal, row, c );

		uint64_t** wedge = wedge_product(gal, b, c, row);

		cout << "wedge product davon ist: " << endl;
		print_matrix(gal, row, row, wedge);
		scalar_matrix(gal, wedge, random_values[i/2], row, row);

		add_matrix(gal, result_matrix, wedge, row, row);
		cout << endl << "Zwischenergebnis für i = " << i << endl;
		print_matrix(gal, row, row, result_matrix);

		//random_value = 1;
		//Result = Result + (*random_values)(i/2) * (M.col(i) * M.col(i+1).t() - M.col(i+1)* M.col(i).t() );

	}

	return result_matrix;
}

int* simple_parity(Galois gal, uint64_t** M, int row, int col){
	Gauss gauss;
	uint64_t* random_values = new uint64_t[col/2]; 		//stores the randomvalues x_i, created in create_Y (needed later on) 
	for (int i = 0; i < col/2; i++){
		random_values[i] = 0;
	}
	int counter = 0;			   						//for positioning vector_parity 

	cout<<"Algorithm start"<< endl;


	uint64_t** Y = create_Y( gal, M, row, col, random_values);

	cout << "Startmatrix M :"<<endl;
	print_matrix(gal,row, col, M);
	cout << "Computed Y :" << endl;
	print_matrix(gal, row, row, Y);

	uint64_t det = gauss.determinant(gal, row, Y);	
	//int rank_Y = rank(Y);
	
	if (det == 0) {
		cout << endl << endl << "THERE IS NO PARITY BASIS" << endl << "Searching for max rank submatrix of Y" << endl;





		//STILL TO BE IMPLEMENTED (Section 6.5 in paper)
		








		/*
		Col<uword> deleted_cols = find_max_submatrix(Y);   //all columns that can be deleted (and still remain full-rank)
	
		deleted_cols = all_except(deleted_cols, Y.n_cols); //all columns that can stay (and still remain full-rank)
		
		M = M.rows(deleted_cols);

		M.print("Reduced M (and new instance for parity-basis-algo with less rows):");
		*/
	}

	Y = create_Y( gal, M, row, col, random_values);

	int* parity_basis = new int[row];		    //in here we will store the indices of the columns that built the parity basis
	for (int i = 0; i < row; i++){
		parity_basis[i] = 0;
	}

	for(int i = 0; i < col; i = i+2){

		//mat Y_prime = Y - random_values(i/2) *  (M.col(i) * M.col(i+1).t() - M.col(i+1)* M.col(i).t() );



		uint64_t** Y_prime= new uint64_t*[row];
		for (int i=0; i<row; i++){
    		Y_prime[i]= new uint64_t[row];
		}
 	
		
		
		for(int row_index = 0 ; row_index < row; row_index++){		//row ist die dimension von Y (Y is row x row matrix, want to update Y)
			for (int col_index = 0; col_index < row ; col_index++){
				//want to add x_i (M[][i] ^ M[][i+1]) to Y
				Y_prime[row_index][col_index] =gal.add(Y[row_index][col_index], gal.multiply(random_values[i/2], gal.add(gal.multiply(M[row_index][i], M[col_index][i+1]), gal.multiply(M[row_index][i+1], M[col_index][i] ) )  ));
			}
		}

		uint64_t det_y_prime = gauss.determinant(gal, row, Y_prime);
		//jetzt ist y_prime aber kaputt, neu berechnen

		for(int row_index = 0 ; row_index < row; row_index++){		//row ist die dimension von Y (Y is row x row matrix, want to update Y)
			for (int col_index = 0; col_index < row ; col_index++){
				//want to add x_i (M[][i] ^ M[][i+1]) to Y
				Y_prime[row_index][col_index] =gal.add(Y[row_index][col_index], gal.multiply(random_values[i/2], gal.add(gal.multiply(M[row_index][i], M[col_index][i+1]), gal.multiply(M[row_index][i+1], M[col_index][i] ) )  ));
			}
		}


		if (det_y_prime != 0){ //if det(Y_prime) != 0
			
			Y = Y_prime;
		}
		else{
			
			parity_basis[counter] = i;
			parity_basis[counter+1] = i+1;
			counter = counter + 2;
		}
	}
	return parity_basis;
}


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

	return add_matrix(gal,bct, cbt, size, size);

 	

}


/*
 * computes mat = scalar * mat
 */
uint64_t** scalar_matrix(Galois gal, uint64_t** mat, uint64_t scalar, int row, int col){
	for(int i = 0; i < col ; i++){
		for(int j = 0; j < row; j++){
			mat[i][j] = gal.multiply(mat[i][j], scalar);
		}
	}
	return mat;
}

/*
 * computes mat1 = mat1 + mat2
 */
uint64_t** add_matrix(Galois gal, uint64_t** mat1, uint64_t** mat2, int row, int col){
	for(int i = 0; i < col ; i++){
		for(int j = 0; j < row; j++){
			mat1[i][j] = gal.add(mat1[i][j], mat2[i][j]);
		}
	}
	return mat1;
}



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

void print_vector(Galois g, int size, uint64_t* vector)
{
  cout<<"Vector is:" <<endl;
  for (int i = 0; i < size; i++){
  	
  		cout << g.to_string(vector[i]) << endl;
 
  }
  
}


void print_vector_normal(int size, int* vector)
{
  cout<<"Vector is:" <<endl;
  for (int i = 0; i < size; i++){
  	
  		cout <<(vector[i]) << " ";
   
  }
  cout <<endl;
  
}



int main(){
	
	Galois gal;
  	Gauss gauss;
  	gal.set_w(8);
  	gal.set_mode_naive();
  	gal.seed();
/*
  	int row = 4;
  	int col = 6;

  	uint64_t[row][row] matrix**;

  	for(int i = 0; i < row; i++){
  		for(int j = 0 ; j < row; j++){
  			matrix[i][j] = gal.uniform_random_element();
  		}
  	}
  	gauss.determinant(gal, row, matrix);

  	uint64_t zahl1 = gal.uniform_random_element();
  	uint64_t zahl2 = gal.uniform_random_element();
  	
	return 0;
*/

	
	int row = 4;
	int col = 6;
	
 
 	//construct a 4x6 matrix
	uint64_t** matrix= new uint64_t*[row];
 
	for (int i=0; i<row; i++){
    	matrix[i]= new uint64_t[col];
	}
 	
	//fill it with random numbers
	for (int i=0; i<row; i++){
    	for(int j=0; j<col; j++){
        	
        	matrix[i][j]= gal.uniform_random_element();
        	
 			
    	}
	}



/*

	matrix[0][0] = (uint64_t) 0;
	matrix[0][1] = (uint64_t) 3;
	matrix[0][2] = (uint64_t) 7;
	matrix[0][3] = (uint64_t) 6;
	matrix[0][4] = (uint64_t) 13;
	matrix[0][5] = (uint64_t) 2;

	matrix[1][0] = (uint64_t) 11;
	matrix[1][1] = (uint64_t) 12;
	matrix[1][2] = (uint64_t) 13;
	matrix[1][3] = (uint64_t) 1;
	matrix[1][4] = (uint64_t) 0;
	matrix[1][5] = (uint64_t) 7;

	matrix[2][0] = (uint64_t) 14;
	matrix[2][1] = (uint64_t) 15;
	matrix[2][2] = (uint64_t) 0;
	matrix[2][3] = (uint64_t) 12;
	matrix[2][4] = (uint64_t) 1;
	matrix[2][5] = (uint64_t) 4;

	matrix[3][0] = (uint64_t) 5;
	matrix[3][1] = (uint64_t) 6;
	matrix[3][2] = (uint64_t) 4;
	matrix[3][3] = (uint64_t) 12;
	matrix[3][4] = (uint64_t) 11;
	matrix[3][5] = (uint64_t) 5;

	matrix[4][0] = (uint64_t) 11;
	matrix[4][1] = (uint64_t) 3;
	matrix[4][2] = (uint64_t) 8;
	matrix[4][3] = (uint64_t) 6;
	matrix[4][4] = (uint64_t) 1;
	matrix[4][5] = (uint64_t) 3;

	matrix[5][0] = (uint64_t) 5;
	matrix[5][1] = (uint64_t) 13;
	matrix[5][2] = (uint64_t) 14;
	matrix[5][3] = (uint64_t) 10;
	matrix[5][4] = (uint64_t) 12;
	matrix[5][5] = (uint64_t) 3;
*/
	//int help = mat_lin_dep(gal, matrix, row);
	//cout << "help ist" << help;
	


	print_matrix(gal,row, col, matrix);
		
	//uint64_t test = gauss.determinant(gal, row, matrix);

	//cout << "det ist hier " << gal.to_string(test) << endl << endl;

/*
	uint64_t* vector = new uint64_t[col/2];
	cout << "row ist " << row << "und col ist " <<col << endl;
	uint64_t** Y = create_Y( gal, matrix, row, col, vector);
	cout << "Y = " << endl;
	print_matrix(gal, row, row , Y);

	uint64_t det = gauss.determinant(gal, row, Y);

	cout << "Det ist " << gal.to_string(det);
*/	
	
	int* result = simple_parity(gal, matrix, row, col);

	cout << "ENDERGEBNIS " << endl;
	print_vector_normal(row, result);


/*
	uint64_t zahl1 = gal.uniform_random_element();
  	uint64_t zahl2 = gal.uniform_random_element();

  	cout << endl << endl << gal.to_string(zahl1) << " mal " << gal.to_string(zahl2) << " ist " << gal.to_string(gal.multiply(zahl1,zahl2));
*/

  	//int help = mat_lin_dep(gal, Y, row);
  	//cout << "help ist" << help;
  
	return 0;

}