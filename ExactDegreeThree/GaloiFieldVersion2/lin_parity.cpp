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
mat create_Y(mat M, vec *random_values, int max_random_value){
	mat Result(M.n_rows, M.n_rows, fill::zeros);

	for(int i = 0 ; i < M.n_cols ; i = i+2){

		if ( (*random_values)(i/2) == 0 ) (*random_values)(i/2) = get_random_value(max_random_value); // if no specific random_value vector is handed, a new random value is assigned 	
		cout << "Random value ist " << (*random_values)(i/2) << endl;									//if something goes wrong, you can see the random value on the console
		//random_value = 1;
		Result = Result + (*random_values)(i/2) * (M.col(i) * M.col(i+1).t() - M.col(i+1)* M.col(i).t() );

	}

	return Result;
}


/*
 * returns a size-1 dimensional vector (0,1,...,i-1,i+1,..size)
 */
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
Col<uword> find_max_submatrix(mat Y){
	int counter = 0; //counter for the amount of deleted columns
	int rank_Y = matRank(Y);

	Col<uword> deleted_cols(Y.n_cols - rank_Y);
	int pos = 0;
	int columns = Y.n_cols;
	for(int i = 0; i < columns; i++){
		
		mat help = delete_column(i-counter, Y);
	
		if (matRank(help) == rank_Y){
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
std::vector<int> simple_parity(mat M, int max_random_value){
	vec random_values(M.n_cols/2, fill::zeros); 		//stores the randomvalues x_i, created in create_Y (needed later on) 
	int counter = 0;			   						//for positioning vector_parity 

	cout<<"Algorithm start"<< endl;
	mat Y = create_Y(M, &random_values, max_random_value);

	M.print("Startmatrix M :");
	Y.print("Computed Y :");

	int rank_Y = rankMat(Y);
	
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
