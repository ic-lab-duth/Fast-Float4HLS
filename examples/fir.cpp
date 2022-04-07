#include "fast_float.h"
#include "mc_scverify.h"

// Implementation of an FIR filter using the fast-float library
#pragma hls_design
template<int N, typename T = ffp32>
void CCS_BLOCK(fir)(T inp, T coeff[N], T &out) {
    
    static T x[N];
    T sum;

    #pragma hls_unroll
    for (int i=N-1; i>0; i--) {
      x[i] = x[i-1];  
    }
    x[0] = inp;

    sum = 0;
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

CCS_MAIN(int, char **) {
    
    typedef ffp32 DATA_T; 
     
    const int TAPS = 5;
    const int ITER = 10;

    float refIn, refCoeff[TAPS], refOut;
    DATA_T In, Coeff[TAPS], Out;
    
    float HI = 65415456364.0;
    float LO = -654364.0;

    srand(time(NULL));

    for (int i=0; i<TAPS; i++) {
      refCoeff[i] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO))); 
      Coeff[i] = refCoeff[i];
    }

    for (int i=0; i<ITER; i++) {
      refIn =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
      In = refIn;

      fir<TAPS,DATA_T>(In,Coeff,Out);
      refFir<TAPS>(refIn,refCoeff,refOut);

      std::cout << "OUT: " << Out.to_float() << std::endl;
      std::cout << "REF: " << refOut << "\n" << std::endl;
    }
        
    CCS_RETURN(0);
}
