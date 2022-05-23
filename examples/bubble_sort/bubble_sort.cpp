// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Bubble Sort sorting algorithm
#include <iostream>
#include <random>
#include "fast_float.h"
#include "mc_scverify.h"

typedef ffp32 T;

#pragma hls_design
template<int N>
void CCS_BLOCK(bubbleSort)(T inA[N]) {
    
  for (int i=0; i<N-1; i++){
    for (int j=0; j<N-i-1; j++) {
      if (inA[j] > inA[j+1]) {
        T temp = inA[j];
        inA[j] = inA[j+1];
        inA[j+1] = temp;
      } 
    }
  }
};

float random_float() {
  float HI = rand();
  float LO = -rand();
  return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
};

CCS_MAIN(int, char**) {
  
  const int NSIZE = 10;
  T inA[NSIZE];

  srand(time(NULL));

  for (int i=0; i<NSIZE; i++) {
    inA[i] =  random_float();
    std::cout << inA[i].to_float() << std::endl;
  }
  std::cout << std::endl;
  
  bubbleSort<NSIZE>(inA);

  for (int i=0; i<NSIZE; i++) {
    std::cout << inA[i].to_float() << std::endl;
  }

  CCS_RETURN(0);
}
