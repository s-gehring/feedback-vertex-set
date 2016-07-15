#pragma once
#include <string>
#include <random>

class Galois
{
  public:
    std::mt19937 gen;
    int w;
    int mode;
    uint64_t add(uint64_t, uint64_t);
    uint64_t multiply(uint64_t, uint64_t);
    uint64_t naive_multiply(uint64_t, uint64_t);
    uint64_t table_multiply(uint64_t, uint64_t);
    uint64_t logtb_multiply(uint64_t, uint64_t);
    uint64_t clmul_multiply(uint64_t, uint64_t);
    uint64_t divide(uint64_t, uint64_t);
    uint64_t naive_divide(uint64_t, uint64_t);
    uint64_t table_divide(uint64_t, uint64_t);
    uint64_t logtb_divide(uint64_t, uint64_t);
    uint64_t clmul_divide(uint64_t, uint64_t);
    uint64_t inverse(uint64_t);
    uint64_t log(uint64_t);
    uint64_t ilog(uint64_t);
    uint64_t uniform_random_element();
    void seed();
	void seed(int se);

    void set_mode_naive();
    void set_mode_table();
    void set_mode_logtb();
    void set_mode_pcmul();
    void create_table();
    void create_logtb();
    void set_w(int);
    std::string to_string(uint64_t);

    static const int NAIVE = 0;
    static const int TABLE = 1;
    static const int LOGTB = 2;
    static const int CLMUL = 3;

    uint64_t** m_table;
    uint64_t** d_table;
    uint64_t*  logtb;
    uint64_t*  ilogtb;
};
