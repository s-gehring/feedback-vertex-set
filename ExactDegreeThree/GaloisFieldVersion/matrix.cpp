#include "galois.h"
#include "gauss.h"
#include "lin_parity.h"
#include "matrix.h"


/*
 * computes mat = scalar * mat
 */
void scalar_matrix(Galois gal, uint64_t** mat, uint64_t scalar, int row, int col){
	for(int i = 0; i < col ; i++){
		for(int j = 0; j < row; j++){
			mat[i][j] = gal.multiply(mat[i][j], scalar);
		}
	}

}

/*
 * computes mat1 = mat1 + mat2
 */
void add_matrix(Galois gal, uint64_t** mat1, uint64_t** mat2, int row, int col){
	for(int i = 0; i < col ; i++){
		for(int j = 0; j < row; j++){
			mat1[i][j] = gal.add(mat1[i][j], mat2[i][j]);
		}
	}
	
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


/*
    Invertiert eine NxN Matrix mit Hilfe des Gauß-Jordan-Algorithmus.
    mat - Matrix die Invertiert werden soll.
    inv - Die Inverse der Matrix mat.
    return: invertierte matrix
*/
uint64_t** invertMatrix(Galois gal, uint64_t** mat, int size)
{

	uint64_t** inv= new uint64_t*[size];
	for (int i=0; i<size; i++){
    	inv[i]= new uint64_t[size];
	}
    // Eine sizex2*size Matrix für den Gauß-Jordan-Algorithmus aufbauen
   
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
 
    // Gauß-Algorithmus.
    for(int k = 0; k < size-1; ++k)
    {
        // Zeilen vertauschen, falls das Pivotelement eine Null ist
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
                    return 0; // Es gibt kein Element != 0
            }
        }
 
        // Einträge unter dem Pivotelement eliminieren
        for(int i = k+1; i < size; ++i)
        {
            uint64_t p =  gal.divide(A[i][k],A[k][k]);     
            for(int j = k; j < 2*size; ++j)
           
            	A[i][j] = gal.add(A[i][j], gal.multiply(A[k][j], p) );
        }
    }
 
    // Determinante der Matrix berechnen
   
    uint64_t det = A[0][0];
  	for (int i = 1; i < size; i++)
  	{
    det = gal.multiply(det, A[i][i]);
  	}


    if(det == 0)  // Determinante ist =0 -> Matrix nicht invertierbar
        return 0;
 
    // Jordan-Teil des Algorithmus durchführen
    for(int k = size-1; k > 0; --k)
    {
        for(int i = k-1; i >= 0; --i)
        {
            uint64_t p = gal.divide(A[i][k],A[k][k]);
            for(int j = k; j < 2*size; ++j)
                
            	A[i][j] = gal.add(A[i][j], gal.multiply(A[k][j], p));
        }
    }
 
    // Einträge in der linker Matrix auf 1 normieren und in inv schreiben
    for(int i = 0; i < size; ++i)
    {
        uint64_t f = A[i][i];
        for(int j = size; j < 2*size; ++j)
            inv[i][j-size] = gal.divide(A[i][j],f);
    }
 	free(A);
    return inv;
}

/*
 * Matrix multiplication
 * creates a new matrix C = A * B (so col_A should be row_B, for clarity bith are needed)
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


/*
 * Berechnet I + V^T M U mit V 2xn und U nx2-matrizen und M n x n Matrix
 *Das wird in der small rank update formula von Sherman-Morrison-Woddbury benutzt
 */
uint64_t** SMW_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** U, int size){

	

	uint64_t** VM_hilf = multiplication_matrix(gal, V, 2, size, M, size, size);

	uint64_t** result = multiplication_matrix(gal, VM_hilf, 2, size, U, size, 2);

	free(VM_hilf);

	//now add 1 to the diagonal
	result[0][0] = gal.add( result[0][0], (uint64_t) 1);
	result[1][1] = gal.add( result[0][0], (uint64_t) 1);
	return result;
}

/*
 * Berechnet M U (I + V M U) V M   (M = M^¹ und V = V_transponiert 2xn)
 * Benötigt um M⁻¹ upzudaten 
 * NEEDED FOR SMALL RANK UPDATE 4.1
 */
uint64_t** SMW_inverse_update_matrix(Galois gal, uint64_t** V, uint64_t** M, uint64_t** U, int size ){

	uint64_t** MU = multiplication_matrix(gal, M, size, size, U, size, 2); //size x 2
	uint64_t** SMW = SMW_matrix(gal, V, M, U, size);					   // 2 x 2
	uint64_t** VM = multiplication_matrix(gal, V, 2, size, M, size, size); // 2 x size

	uint64_t** MU_SMW = multiplication_matrix(gal, MU, size, 2, SMW, 2, 2); //size x 2
	uint64_t** result = multiplication_matrix(gal, MU_SMW, size, 2, VM, 2, size); //size x size

	free (MU);
	free(SMW);
	free(VM);
	free(MU_SMW);

	return result;
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
	add_matrix(gal,bct, cbt, size, size);
	free(cbt);
	return bct;

}