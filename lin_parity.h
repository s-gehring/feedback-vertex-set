#ifndef LIN_PARITY_H_
#define LIN_PARITY_H_

#include <iostream>
#ifndef debug
    #define debug if(1)
#endif
#include <cstdlib> 
#include <ctime>
#include "galois.h"


using namespace std;

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
uint64_t** create_Y(Galois & gal, uint64_t **M, int row, int col, uint64_t *random_values);


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
uint64_t** get_U(Galois & gal, uint64_t** M, int size, uint64_t random, int i);


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
 uint64_t** get_V(Galois & gal, uint64_t** M, int size, int i);



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
int* simple_parity_fast(Galois & gal, uint64_t** M, int row, int col, int* length);




#endif /* LIN_PARITY_H_ */
