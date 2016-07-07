#ifndef LIN_PARITY_H_
#define LIN_PARITY_H_

#include <iostream>

#include <cstdlib> 
#include <ctime>
#include "galois.h"



//using namespace arma;

using namespace std;
//using namespace boost::numeric::ublas;
int get_random_value(int max);
/* Creates the compact Matrix Y and stores the random values x_i in vector random_values
 * so random_values has to be a "column of M-half-dim" vec (a pointer of it)
 * random values are integral between 1 and max_random_value
 */ 
uint64_t** create_Y(Galois gal, uint64_t **M, int row, int col, uint64_t *random_values);

/*
 * basically not needed, this is a computation with 2 intermediate steps (but this means there a two row x row matrices)
 * in the other computation (ref create_Y() ) there is a computation for every cell in one line
 * For debugging very helpful
 */
uint64_t** create_Y_naive(Galois gal, uint64_t **M, int row, int col, uint64_t *random_values);


/*
 * This implementations does not use the Sherman Morrisin Woodbury fast rank update formula
 */
int* simple_parity(Galois gal, uint64_t** M, int row, int col);

/*
 * gibt U = x_i * ( b_i  c_i ) zurück, also x_i * ( i-th col   i+1-th col )
 * U ist also size x 2 matrix
 */
uint64_t** get_U(Galois gal, uint64_t** M, int size, uint64_t random, int i);

/*
 * gibt V = ( -c_i  b_i )^T zurück, also ( - i+1-th col   i-th col) ^T
 * V ist also 2 x size matrix
 */

 uint64_t** get_V(Galois gal, uint64_t** M, int size, int i);



/*
 * This implemetation uses the small rank update formula and also return a maximum number of pairs, if there is no
 * parity basis
 * in length it returns the number of lin. independent vectors (so there are lenght/2 pairs)
 */
int* simple_parity_fast(Galois gal, uint64_t** M, int row, int col, int* length);
//std::vector<int> simple_parity(mat M, int max_random_value);



#endif /* LIN_PARITY_H_ */
