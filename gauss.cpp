#include <iostream>
#include "gauss.h"

using namespace std;
 /**
	* @brief calculates the determinat
	*/
uint64_t Gauss::determinant(Galois g, int size, uint64_t** matrix)
{
  uint64_t** transform = upper_triangle_transform(g, size, matrix);
  uint64_t result = transform[0][0];
  for (int i = 1; i < size; i++)
  {
    result = g.multiply(result, transform[i][i]);
  }
  return result;
}

 /**
	* @brief transforms the matrix to upper triangle form
	*/
uint64_t** Gauss::upper_triangle_transform(Galois g, int size, uint64_t** matrix)
{

  for (int n = 0; n < size; n++)
  {
    int nonzero = -1;
    bool found = false;
    int k = n;

    while (!found && k < size)
    {
      if (matrix[k][n] != 0)
      {
        nonzero = k;
	found = true;
      }
      k++;
    }
     
    if  (found) //or column is 0
    {  
      if (nonzero != n)
      {

        for (int j = n; j < size; j++)
        {
          uint64_t old_n = matrix[n][j];
          uint64_t old_nonzero = matrix[nonzero][j];
          matrix[n][j] = old_nonzero;
          matrix[nonzero][j] = old_n;
        }
  
      }

      for (int i = n+1; i < size; i++)
      {
   
        if (matrix[i][n] != 0)
        {
          uint64_t div = g.divide(matrix[i][n], matrix[n][n]);

          for (int j = n; j < size; j++)
          {
            uint64_t mul = g.multiply(div, matrix[n][j]);
            uint64_t add = g.add(mul, matrix[i][j]);
            matrix[i][j] = add;
          }

        }
      }     
    }
  }
  
  return matrix;
}

