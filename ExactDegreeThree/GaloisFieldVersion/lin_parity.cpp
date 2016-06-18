
#include "galois.h"
#include "gauss.h"
#include "lin_parity.h"
#include "matrix.h"



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
		if (col < 30) cout << "Random value ist " << gal.to_string(random_values[i/2]) << endl;	//if something goes wrong, you can see the random value on the console

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
		//random_values[i/2] = (uint64_t) 1;
		cout << "Random value ist " << gal.to_string(random_values[i/2]) << endl;	//if something goes wrong, you can see the random value on the console

		uint64_t* b = new uint64_t[row];
		
		for (int j = 0; j < row; j++){
			b[j] = M[j][i];
		}

		uint64_t* c = new uint64_t[row];
		for (int j = 0; j < row; j++){
			c[j] = M[j][i+1];
		}

		//cout << "b ist " << endl;
		//print_vector(gal, row, b );
		//cout << "c ist " << endl;
		//print_vector(gal, row, c );

		uint64_t** wedge = wedge_product(gal, b, c, row);

		//cout << "wedge product davon ist: " << endl;
		//print_matrix(gal, row, row, wedge);
		scalar_matrix(gal, wedge, random_values[i/2], row, row);

		cout << "komponente von i = " << i  << endl << endl << endl;
		print_matrix(gal, row, row, wedge);

		add_matrix(gal, result_matrix, wedge, row, row);
		//cout << endl << "Zwischenergebnis für i = " << i << endl;
		//print_matrix(gal, row, row, result_matrix);

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


/*
 * gibt U = x_i * ( b_i  c_i ) zurück, also x_i * ( i-th col   i+1-th col )
 * U ist also size x 2 matrix
 */
uint64_t** get_U(Galois gal, uint64_t** M, int size, uint64_t random, int i){
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

/*
 * gibt V = ( -c_i  b_i )^T zurück, also ( - i+1-th col   i-th col) ^T
 * V ist also 2 x size matrix
 */

 uint64_t** get_V(Galois gal, uint64_t** M, int size, int i){
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

int* simple_parity_fast(Galois gal, uint64_t** M, int row, int col){
	Gauss gauss;
	uint64_t* random_values = new uint64_t[col/2]; 		//stores the randomvalues x_i, created in create_Y (needed later on) 
	for (int i = 0; i < col/2; i++){
		random_values[i] = 0;
	}
	int counter = 0;			   						//for positioning vector_parity 

	cout<<"------------------Algorithm start------------------"<< endl;
	cout << "Startmatrix M :"<<endl;
	if (row < 100) print_matrix(gal,row, col, M);
	else cout << "too big, will not write it down" <<endl;

	cout << "Computing Y...";
	uint64_t** Y = create_Y( gal, M, row, col, random_values);
	cout << "done." << endl;
	//print_matrix(gal, row, row, Y);

	cout << "Computing Inverse Y^-1 ...";
	uint64_t** Y_inverse = invertMatrix(gal, Y,  row);
	cout << "done" << endl;
	

	cout << "Computing determinant of Y ...";
	uint64_t det = gauss.determinant(gal, row, Y);	
	cout << "done. Det(Y) = " << gal.to_string(det) << endl;

	
	
	if (det == 0) {
		//Hier ist Y zerstört durch det
		cout << endl << endl << "THERE IS NO PARITY BASIS" << endl << "Searching for max rank submatrix of Y" << endl;
		for(int k = 0; k < 1000000; k++){
			cout << "ABORT" << endl;
		}
		/*
		Col<uword> deleted_cols = find_max_submatrix(Y);   //all columns that can be deleted (and still remain full-rank)
	
		deleted_cols = all_except(deleted_cols, Y.n_cols); //all columns that can stay (and still remain full-rank)
		
		M = M.rows(deleted_cols);

		M.print("Reduced M (and new instance for parity-basis-algo with less rows):");
		*/
	}
	cout << "Computing again Y...";
	Y = create_Y( gal, M, row, col, random_values);
	cout << "done." << endl;
	//print_matrix(gal, row, row, Y);



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

		
		uint64_t** V = get_V(gal, M, row, i);
		uint64_t** U = get_U(gal, M, row,random_values[i/2], i);

		
 		
 		uint64_t** SMW_det = SMW_matrix(gal,  V, Y_inverse, U, row);
		uint64_t** SMW = copy_matrix(SMW_det, 2,2);
		uint64_t det_SMW = gauss.determinant(gal, 2, SMW_det);

		free(SMW_det);

		cout<< "det von Y_prime durch trick ist " << gal.to_string(det_SMW) << endl;
		


		if (det_SMW != 0){ //if det(Y_prime) != 0

			//Setzt Y auf Y' ( also Y + UV)
			uint64_t** hilf = multiplication_matrix(gal, U, row, 2, V, 2, row);
			add_matrix(gal, Y, hilf, row, row);

			//update Y_inverse
			uint64_t** update_matrix = SMW_inverse_update_matrix( gal,  V,  Y_inverse,  U, row );
			add_matrix( gal, Y_inverse, update_matrix,  row, row);
		
		}
		else{
			
			parity_basis[counter] = i;
			parity_basis[counter+1] = i+1;
			counter = counter + 2;
		}
	}
	return parity_basis;
}




/*
 * returns a size-1 dimensional vector (0,1,...,i-1,i+1,..size)
 */
 /*
Col<uword> every_column_except(int col, int size){//vector muss <uword> sein, um in funktion Y.cols(hilfsvec) eingesetzt werden zu können
	Col<uword> vector(size-1);
	int counter = 0;
	for (int i = 0; i < size; i++ ){
		if ( i != col){
		  	vector(counter) = i;
		  	counter++;
		}
	}
	return vector;
}


mat delete_column(int col, mat matrix){
	Col<uword> vector = every_column_except(col, matrix.n_cols);
	return matrix.cols(vector);
}


/*
 * return the vector of indices, which, if we delete all the columns indices with the return indices, is still of full rank
 * therefor it goes through all columns, checks if we can delete this column without decreasing the rank
 */
 /*
Col<uword> find_max_submatrix(mat Y){
	int counter = 0; //counter for the amount of deleted columns
	int rank_Y = arma::rank(Y);

	Col<uword> deleted_cols(Y.n_cols - rank_Y);
	int pos = 0;
	int columns = Y.n_cols;
	for(int i = 0; i < columns; i++){
		
		mat help = delete_column(i-counter, Y);
	
		if (arma::rank(help) == rank_Y){
			//we can delete column i without loosing rank
			//store this information
			deleted_cols(pos) = i;  
			pos++;
			//update the Matrix
			Y = help;
			counter++;  //store the information, that we deleted a column (so the next col to check is (i-counter))
		}

		//otherwise if the rank differs, we do nothing and check the next column

	}
	return deleted_cols;

}


/*
 * input is a vector of indices 
 * outputs a vector of all indices ascending from 0 to size-1 without those in vec_not
 */
 /*
Col<uword> all_except(Col<uword> vec_not, int size){
	Col<uword> vec_result(size-vec_not.n_elem);
	int pos = 0;			//to count up the vec_not vector (the forbidden indices in vec_not are ascending)
	int counter = 0;		//to count up the vec_result vector

	for (int i = 0; i < size ; i++){
		
		if ((vec_not.n_elem > pos) && (vec_not(pos) == i)){   //the first condition takes care of the problem, that after the last forbidden index, pos is still increased
			//i is skipped, increase pos to the next forbidden index
			pos++;
		}
		else{
			//index i is not forbidden, write it to the result-vertex
			vec_result(counter) = i;
			counter++;
		}
	}
	return vec_result;
}
/*
 * The main algorithm 4.1 with handeliing of no parity-basis (section 6.5)
 */

 /*
std::vector<int> simple_parity(mat M, int max_random_value){
	vec random_values(M.n_cols/2, fill::zeros); 		//stores the randomvalues x_i, created in create_Y (needed later on) 
	int counter = 0;			   						//for positioning vector_parity 

	cout<<"Algorithm start"<< endl;
	mat Y = create_Y(M, &random_values, max_random_value);

	M.print("Startmatrix M :");
	Y.print("Computed Y :");

	int rank_Y = arma::rank(Y);
	
	if (rank_Y != Y.n_cols) {
		cout << endl << endl << "THERE IS NO PARITY BASIS" << endl << "Searching for max rank submatrix of Y" << endl;
		
		Col<uword> deleted_cols = find_max_submatrix(Y);   //all columns that can be deleted (and still remain full-rank)
	
		deleted_cols = all_except(deleted_cols, Y.n_cols); //all columns that can stay (and still remain full-rank)
		
		M = M.rows(deleted_cols);

		M.print("Reduced M (and new instance for parity-basis-algo with less rows):");
	}

	Y = create_Y(M, &random_values, max_random_value);

	//vec parity_basis(M.n_rows, fill::zeros);		    //in here we will store the indices of the columns that built the parity basis
	std::vector<int>parity_basis(M.n_rows,0);

	for(int i = 0; i < M.n_cols; i = i+2){

		mat Y_prime = Y - random_values(i/2) *  (M.col(i) * M.col(i+1).t() - M.col(i+1)* M.col(i).t() );

		if (arma::rank(Y_prime) == Y_prime.n_cols){ //if det(Y_prime) != 0
			Y = Y_prime;
		}
		else{
			//parity_basis(counter) = i;
			//parity_basis(counter+1) = i+1;
			parity_basis[counter]=i;
			parity_basis[counter+1]=i+1;
			counter = counter + 2;
		}
	}
	return parity_basis;
}
*/












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

	
	int row = 400;
	int col = 600;
	
 
 	//construct a row x col matrix
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
	cout<< "Startmatrix = ";
	print_matrix(gal, row, col, matrix);

	/*
	uint64_t* random_values = new uint64_t[col/2]; 		//stores the randomvalues x_i, created in create_Y (needed later on) 
	for (int i = 0; i < col/2; i++){
		random_values[i] = 0;
	}

	uint64_t** Y = create_Y( gal, matrix, row, col, random_values);
	
	//print_matrix(gal, row, row, Y);

	uint64_t** Y_inverse = invertMatrix(gal, Y,  row);
	
	cout << "Y is" << endl;
	print_matrix(gal, row, row, Y);
	cout << "Y-inverse is" <<endl;
	print_matrix(gal, row, row,Y_inverse);

	uint64_t** test = multiplication_matrix(gal, Y, row, row, Y_inverse, row, row);
	cout << "test mult = " << endl;
	print_matrix(gal, row, row, test);


	uint64_t** V = get_V(gal, matrix, row, 0);
	uint64_t** U = get_U(gal, matrix, row,random_values[0/2], 0);

	uint64_t** UV = multiplication_matrix(gal, U, row, 2, V, 2, row);
	cout << "UV =" << endl;
	print_matrix(gal, row, row, UV);
	add_matrix(gal, Y, UV, row, row);

	cout << "new Y = " << endl;
	print_matrix(gal, row, row, Y);



	uint64_t** VY_I = multiplication_matrix(gal, V, 2, row, Y_inverse, row, row);  //2 x row
	uint64_t** VY_IU = multiplication_matrix(gal, VY_I, 2, row, U, row, 2);      //2 x 2
	print_matrix(gal, 2, 2, VY_IU);
	VY_IU[0][0] = gal.add(VY_IU[0][0] , (uint64_t) 1);
	VY_IU[1][1] = gal.add(VY_IU[1][1] , (uint64_t) 1);
	
	print_matrix(gal, 2, 2, VY_IU);

	*/
	
/*
	uint64_t** inv = invertMatrix(gal, matrix,  row);
	//uint64_t det = gauss.determinant(gal, row, matrix);
	print_matrix(gal, row, row, inv);

	uint64_t** matrix3 = multiplication_matrix(gal, matrix, row, row, inv , row, row);

	print_matrix(gal, row, row, matrix3);
*/

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
	



	int* result = simple_parity_fast(gal, matrix, row, col);

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
