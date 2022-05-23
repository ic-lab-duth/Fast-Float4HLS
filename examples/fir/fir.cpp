// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build an FIR filter
#include <iostream>
#include <random>
#include "fast_float.h"
#include "mc_scverify.h"

typedef ffp32 T;

// Implementation of an FIR filter using the fast-float library
#pragma hls_design
#pragma hls_pipeline_init_interval 1
template<int N>
void CCS_BLOCK(fir)(T inp, T coeff[N], T &out) {
    
    static T x[N], sum;

    #pragma hls_unroll
    for (int i=N-1; i>0; i--) {
      x[i] = x[i-1];  
    }
    x[0] = inp;

    sum = 0.0;
    #pragma hls_unroll
    for (int i=0; i<N; i++) {
        x[i].fpma_dual(coeff[i],sum,sum);
    }
    out = sum;
};

// Reference FIR filter using C++ float datatype
template<int N>
void refFir(float inp, float coeff[N], float &out) {
    
    static float x[N];
    float sum;

    for (int i=N-1; i>0; i--) {
      x[i] = x[i-1];  
    }
    x[0] = inp;

    sum = 0;
    for (int i=0; i<N; i++) {
        sum += x[i] * coeff[i];
    }
    out = sum;
};

float random_float() {
  float HI = rand();
  float LO = -rand();
  return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
};

CCS_MAIN(int argc, char * argv[]) {
        
    const int TAPS = 5;
    const int ITER = 10;

    float refIn, refCoeff[TAPS], refOut;
    T In, Coeff[TAPS], Out;
    
    srand(time(NULL));

    for (int i=0; i<TAPS; i++) {
      refCoeff[i] = random_float(); 
      Coeff[i] = refCoeff[i];
    }

    for (int i=0; i<ITER; i++) {
      refIn =  random_float();
      In = refIn;

      fir<TAPS>(In,Coeff,Out);
      refFir<TAPS>(refIn,refCoeff,refOut);

      std::cout << "OUT: " << Out.to_float() << std::endl;
      std::cout << "REF: " << refOut << "\n" << std::endl;
    }
        
    CCS_RETURN(0);
}
