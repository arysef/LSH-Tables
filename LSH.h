#pragma once

#include "Reservoir.h"
#include <iostream>
#include <omp.h>
#include <vector>

class LSH {
  private:
    unsigned int L;
    unsigned int reservoir_size;
    unsigned int range_pow;
    unsigned int range;
    Reservoir **reservoirs;

  public:
    LSH();




    //num_items is the number of items we're inserting
    // items is a pointer to the first item in an array of tems.
    // hashes is a pointer to an array of all the k*L hashes. 
    //puts in all the hashes we need to hash 
    void insert(unsigned int num_items, unsigned int *items, unsigned int *hashes);

    void insert(unsigned int item, unsigned int *hashes);

    void retrieve(unsigned int num_query, unsigned int *hashes, unsigned int *results_buffer);

    void top_k(unsigned int num_query, unsigned int top_k, unsigned int *hashes,
               unsigned int *selection);

    void reset();

    void view();

    void add_random_items(unsigned int num_items, bool verbose);

    ~LSH();
};