#ifndef MATRIX_H_
#define MATRIX_H_


#include <iostream>

#include <cstdlib> 
#include <ctime>


/*
 * computes mat = scalar * mat
 */
void scalar_matrix(Galois gal, uint64_t** mat, uint64_t scalar, int row, int col);

/*
 * computes mat1 = mat1 + mat2
 */
void add_matrix(Galois gal, uint64_t** mat1, uint64_t** mat2, int row, int col);

void print_matrix(Galois g, int size_row, int size_col, uint64_t** matrix);




void print_vector(Galois g, int size, uint64_t* vector);


void print_vector_normal(int size, int* vector);


/*
*
    Vertauscht zwei Zeilen in einer row x col Matrix.
    line1, line2 - Index der Zeilen, die vertauscht werden sollten.
    return: false, falls line1 oder line2 nicht in der Matrix liegen.
*/
bool swapLine(uint64_t** mat, int row, int col, int line1, int line2);
/*
    Invertiert eine NxN Matrix mit Hilfe des Gauß-Jordan-Algorithmus.
 
    mat - Matrix die Invertiert werden soll.
    inv - Die Inverse der Matrix mat.
    return: invertierte Matrix
*/
uint64_t** invertMatrix(Galois gal, uint64_t** mat, int size);


/*
 * Matrix multiplication
 * creates a new matrix C = A * B (so col_A should be row_B, for clarity bith are needed)
 */
uint64_t** multiplication_matrix(Galois gal, uint64_t** A, int row_A, int col_A, uint64_t** B, int row_B, int col_B);

/*
 * Berechnet I + V^T M U mit V 2xn und U nx2-matrizen und M n x n Matrix
 *Das wird in der small rank update formula von Sherman-Morrison-Woddbury benutzt
 * NEEDED FOR SMALL RANK UPDATE 4.1
 */
uint64_t** SMW_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** U, int size);

/*
 * Berechnet M U (I + V M U) V M   (M = M^¹ und V = V_transponiert 2xn)
 * Benötigt um M⁻¹ upzudaten 
  *NEEDED FOR SMALL RANK UPDATE 4.1
 */
uint64_t** SMW_inverse_update_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** U, int size );

uint64_t** wedge_product(Galois gal, uint64_t *b, uint64_t *c, int size);
#endif /* LIN_PARITY_H_ */