
#include "galois.h"
#include "gauss.h"
#include "lin_parity.h"
#include "matrix.h"
#include "matrix2.h"


int get_random_value(int max){
	int random_variable = std::rand();
    return random_variable % max + 1 ;
}


/* Creates the compact Matrix Y and stores the random values x_i in vector random_values
 * so random_values has to be a "column of M-half-dim" vec (a pointer of it)
 * random values are integral between 1 and max_random_value
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
			//cout << "rand ist " << gal.to_string(random_values[i/2]) << endl;
		}
		//random_values[i/2] = (uint64_t) 1;
		//if (col < 30) {debug cout << "Random value ist " << gal.to_string(random_values[i/2]) << endl;}	//if something goes wrong, you can see the random value on the console

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
 * gibt U = x_i * ( b_i  c_i ) zurück, also x_i * ( i-th col   i+1-th col )
 * U ist also size x 2 matrix
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

/*
 * gibt V = ( -c_i  b_i )^T zurück, also ( - i+1-th col   i-th col) ^T
 * V ist also 2 x size matrix
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
	
	//cout << "|";
	//cout.flush();

	
	
	/*cout<<"------------------Algorithm start-----------"<< endl;
	cout << "Startmatrix M :"<<endl;
	if (row < 13) {debug print_matrix(gal,row, col, M);}
	else {cout << "too big, will not write it down" <<endl;}

	cout << "Computing Y...";*/
	uint64_t** Y = create_Y( gal, M, row, col, random_values);

	
	// cout << "done." << endl;
	//print_matrix(gal, row, row, Y);

	uint64_t** Y_copy_det = copy_matrix(Y, row, row);

	// cout << "Computing determinant of Y ...";
	uint64_t det = gauss.determinant(gal, row, Y_copy_det);	
	// cout << "done. Det(Y) = " << gal.to_string(det) << endl;

	my_free(Y_copy_det, row);
	
	if (det == 0) {
		check = 1;
		// cout << endl << endl << "THERE IS NO PARITY BASIS" << endl << "Searching for redundant rows of Y.." ;
		
		std::vector<int> del_row  = findRedundantRows(Y, row, row);
		// cout << "done " <<endl;
		row -= del_row.size(); 
		row_difference = del_row.size();
		sort(del_row.begin(), del_row.end());
		

		//cout <<"row = " << row << "und row_difference = " << row_difference << endl;
		//cout.flush();
		/*for(unsigned int kk = 0; kk < del_row.size(); kk++){
			 cout << "del_row[" << kk << "] = " << del_row[kk] << endl;
		}*/

		if (row == 0){
			//alles wird gelöscht, es gibt kein pair
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


		M = M_prime;


		//cout<<"fedditsch" <<endl;
		//print_matrix(gal, row, col, M);



	}

	


	// cout << "Computing again Y...";
	//my_free(Y, row); 
	if (check == 1){
		my_free(Y, row + row_difference);
		Y = create_Y( gal, M, row, col, random_values);
	}
	// cout << "done." << endl;
	//print_matrix(gal, row, row, Y);

	


	// cout << "Computing Inverse Y^-1 ...";
	uint64_t** Y_inverse = invertMatrix(gal, Y,  row);
	// cout << "done" << endl;

	//print_matrix(gal, row, row, Y_inverse);

	

	parity_basis = new int[row];		    //in here we will store the indices of the columns that built the parity basis
	for (int i = 0; i < row; i++){
		parity_basis[i] = 0;
	}
	
	int pos = 0;			   						//for positioning vector_parity 
	for(int i = 0; i < col; i = i+2){

		//mat Y_prime = Y - random_values(i/2) *  (M.col(i) * M.col(i+1).t() - M.col(i+1)* M.col(i).t() );
		



		uint64_t** V = get_V(gal, M, row, i);
		uint64_t** U = get_U(gal, M, row,random_values[i/2], i);	

					
 	
 		uint64_t** SMW_det = SMW_matrix(gal,  V, Y_inverse, U, row);
 		

		uint64_t** SMW = copy_matrix(SMW_det, 2,2);
		
		uint64_t det_SMW = gauss.determinant(gal, 2, SMW_det);
		
		my_free(SMW_det,2);

		// cout<< "det von Y_prime durch smallrankupdate ist " << gal.to_string(det_SMW) << endl;
			

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
				/*
				cout << "row is " << row << "col is " << col << "pos is " << pos << endl << endl << endl;
				cout << "OVERFLOW STARTE ERNEUT"<<endl <<endl;
				
				
				cout.flush();
				print_vector(gal, (col/2), random_values);
				*/
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

	//print_vector_normal(row, parity_basis);
 }
 //print_vector_normal(row, parity_basis);
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