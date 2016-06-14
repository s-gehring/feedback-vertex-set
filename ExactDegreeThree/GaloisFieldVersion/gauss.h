#ifndef GAUSS_H
#define GAUSS_H


#include "galois.h"

class Gauss
{
  public:
    uint64_t determinant(Galois, int, uint64_t**);
    uint64_t** upper_triangle_transform(Galois, int, uint64_t**);
};
#endif