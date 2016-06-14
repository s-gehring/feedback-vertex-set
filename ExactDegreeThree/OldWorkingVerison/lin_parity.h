#ifndef LIN_PARITY_H_
#define LIN_PARITY_H_

#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cstdlib> 
#include <ctime>

#include <armadillo>

using namespace arma;

using namespace std;
//using namespace boost::numeric::ublas;
int get_random_value(int max);

/* Creates the compact Matrix Y and stores the random values x_i in vector random_values
 * so random_values has to be a "column of M-half-dim" vec (a pointer of it)
 * random values are integral between 1 and max_random_value
 */ 
mat create_Y(mat M, vec *random_values, int max_random_value);

/*
 * returns a size-1 dimensional vector (0,1,...,i-1,i+1,..size)
 */
Col<uword> every_column_except(int col, int size);

mat delete_column(int col, mat matrix);
/*
 * return the vector of indices, which, if we delete all the columns indices with the return indices, is still of full rank
 * therefor it goes through all columns, checks if we can delete this column without decreasing the rank
 */
Col<uword> find_max_submatrix(mat Y);

/*
 * input is a vector of indices 
 * outputs a vector of all indices ascending from 0 to size-1 without those in vec_not
 */
Col<uword> all_except(Col<uword> vec_not, int size);

/*
 * The main algorithm 4.1 with handeliing of no parity-basis (section 6.5)
 */
std::vector<int> simple_parity(mat M, int max_random_value);



#endif /* LIN_PARITY_H_ */
