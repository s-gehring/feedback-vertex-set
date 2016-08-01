#include <bitset>
#include <iostream>
#include "galois.h"
#include <string>
#include <random>
#include "static.cpp"

using namespace std;
 /**
	* @brief add two elements of galois fields
	*/
uint64_t Galois::add(uint64_t arg1, uint64_t arg2)
{
  return arg1 ^ arg2;
}
 /**
	* @brief multiply two elements of galois fields
	*/
uint64_t Galois::multiply(uint64_t arg1, uint64_t arg2)
{
  switch(mode)
  {
    case NAIVE: return naive_multiply(arg1, arg2); break;
    case TABLE: return table_multiply(arg1, arg2); break;
    case LOGTB: return logtb_multiply(arg1, arg2); break;
    case CLMUL: return clmul_multiply(arg1, arg2); break;
    default: return -1;
  }
}

uint64_t Galois::naive_multiply(uint64_t arg1, uint64_t arg2)
{
  uint64_t a = arg1;
  uint64_t b = arg2; 
  uint64_t p = 0;

  for (int i = 0; i < w; i++)
  {
    if (b & 1)
    {
      p ^= a; 
    }
    b >>= 1;

    bool carry = false;
    if (a & (pw[w-1]))
    {
      carry = true;
    }
    a <<= 1;

    if (carry)
    {
      a ^= (primpoly[w]);
    }
  }

  return p;
}

uint64_t Galois::table_multiply(uint64_t arg1, uint64_t arg2)
{
  return m_table[arg1][arg2];
}

uint64_t Galois::logtb_multiply(uint64_t arg1, uint64_t arg2)
{
  if (arg1 == 0 || arg2 == 0) return 0;
  uint64_t sum = logtb[arg1] + logtb[arg2];
  if (sum >= pwm1[w]) sum -= pwm1[w];
  return ilogtb[sum];
}

uint64_t Galois::clmul_multiply(uint64_t arg1, uint64_t arg2)
{
  uint64_t mask = pwm1[32];
  uint64_t x1x0, x2, x3;
  
  __asm__ ( 	//carry-less multiplication
		"movd %3, %%xmm1;" 
	 	"movd %4, %%xmm2;"
		"pclmulqdq $0, %%xmm1, %%xmm2;" 

		//fetch first 64 bits
		"movd %%xmm2, %0;"
                
		//fetch 32 more bits
                "psrldq $8, %%xmm2;"
                "movq %%xmm2, %%xmm3;"
                "movd %5, %%xmm1;"
                "pand %%xmm1, %%xmm3;"
                "movd %%xmm3, %1;"

		//fetch last 32 bits
                "psrldq $4, %%xmm2;"
		"movd %%xmm2, %2;"          

		: "=r" (x1x0), "=r" (x2), "=r" (x3)
		: "r" (arg1), "r" (arg2), "r" (mask)		
	);

  //modulo primitive polynomial
  uint64_t x3d = (x3 << 32) ^ (x3 >> 31) ^ (x3 >> 29) ^ (x3 >> 28) ^ x2;
  uint64_t h1h0 = (x3d << 1) ^ (x3d << 3) ^ (x3d << 4) ^ x3d;
  return h1h0 ^ x1x0;
}

uint64_t Galois::divide(uint64_t arg1, uint64_t arg2)
{
  switch(mode)
  {
    case NAIVE: return naive_divide(arg1, arg2); break;
    case TABLE: return table_divide(arg1, arg2); break;
    case LOGTB: return logtb_divide(arg1, arg2); break;
    case CLMUL: return clmul_divide(arg1, arg2); break;
    default: return -1;
  }

}

uint64_t Galois::naive_divide(uint64_t arg1, uint64_t arg2)
{
  if (arg1 == 0) return 0;
  if (arg2 == 0) return -1;
  return multiply(arg1, inverse(arg2));
}

uint64_t Galois::table_divide(uint64_t arg1, uint64_t arg2)
{
  return d_table[arg1][arg2];
}

uint64_t Galois::logtb_divide(uint64_t arg1, uint64_t arg2)
{
  if (arg1 == 0) return 0; 
  if (arg2 == 0) return -1; 
  uint64_t diff = logtb[arg1]-logtb[arg2];
  if (logtb[arg1] < logtb[arg2]) diff += pwm1[w]; 
  return ilogtb[diff];
}

uint64_t Galois::clmul_divide(uint64_t arg1, uint64_t arg2)
{
  return naive_divide(arg1, arg2);
}

uint64_t Galois::inverse(uint64_t arg1)
{
  //TODO:fix for small w
  uint64_t square = multiply(arg1, arg1); 

  uint64_t result = square;
  for (int i = 0; i < (w-2); i++)
  {
    square = multiply(square, square);
    result = multiply(square, result);   
  }
  return result;

}

uint64_t Galois::log(uint64_t arg1)
{ 
  return logtb[arg1];
}

uint64_t Galois::ilog(uint64_t arg1)
{
  return ilogtb[arg1];
}
 /**
	* @brief set a random seed
	*/
void Galois::seed()
{
  random_device rd;  
  gen.seed(rd());   
}
 /**
	* @brief set a specific seed
	*/
void Galois::seed(int se)
{
	gen.seed(se);
}

 /**
	* @brief generate a random element
	*/
uint64_t Galois::uniform_random_element()
{ 
  int bits = w;
  uint64_t rand = 0;

  uniform_int_distribution<> dist16(0,pwm1[16]);

  while (bits > 16){
    rand <<= 16;
    rand ^= dist16(gen);
    bits -= 16;
  }

  uniform_int_distribution<> distb(0,pwm1[bits]);
  rand <<= bits;
  rand ^= distb(gen);

  return rand;
}

void Galois::set_mode_naive()
{
  mode = NAIVE;
}

void Galois::set_mode_table()
{
  create_table();  
  mode = TABLE;
}

void Galois::create_table()
{
  uint64_t size = pw[w];
  m_table = new uint64_t*[size];
  d_table = new uint64_t*[size];
      
  for (std::size_t i = 0; i < size; i++)
  {
    m_table[i] = new uint64_t[size];
    d_table[i] = new uint64_t[size];
    for (std::size_t j = 0; j < size; j++)
    {
      m_table[i][j] = naive_multiply(i, j);
      d_table[i][j] = naive_divide(i, j);
    }
  }
}

void Galois::set_mode_logtb()
{ 
  create_logtb();
  mode = LOGTB;
}

void Galois::create_logtb()
{
  int size = pw[w];
  logtb = new uint64_t[size];
  ilogtb = new uint64_t[size];

  uint64_t b = 1;

  for (int log = 0; log < (size-1); log++)
  {
    logtb[b] = (unsigned int) log;
    ilogtb[log] = (unsigned int) b;
    b <<= 1;
    if (b & size)
    {
      b ^= primpoly[w];
    }
  }
}

void Galois::set_mode_pcmul()
{
  mode = CLMUL;
  w = 64;
}

 /**
	* @brief set the size of the galois field
	*/
void Galois::set_w(int exponent)
{
  w = exponent;
}

string Galois::to_string(uint64_t arg1)
{
  std::bitset<64> x(arg1);
  string s = x.to_string();
  return s.substr(64-w, w);
}


