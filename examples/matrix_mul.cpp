#include <iostream>
#include "fast_float.h"
#include "mc_scverify.h"

enum OP_TYPE {SEPARATE, FMA, DOT};
typedef ffp32 T;

#pragma hls_design
template<int N, int M, OP_TYPE OP=SEPARATE>
void CCS_BLOCK(matXvec)(T A[N][M], T V[M], T O[N]) {

  for (int i=0; i<N; i++) {
    T sum = 0.0;
    switch (OP) {
    case SEPARATE:
      for (int j=0; j<M; j++) 
        sum+= A[i][j]*V[j];
      break;

    case FMA:
      for (int j=0; j<M; j++) 
        A[i][j].fpma_dual(V[j],sum,sum);
      break;

    case DOT:
        sum.dotProd<M>(A[i],V);
      break;
    
    default:
      for (int j=0; j<M; j++) 
        sum+= A[i][j]*V[j];
      break;
    }
    O[i] = sum;
  }
}

#pragma hls_design
template<int N, int M, int K, OP_TYPE OP=SEPARATE>
void CCS_BLOCK(matXmat) (T A[N][M], T B[M][K], T O[N][K]) {
  for (int i=0; i<K; i++) {
    T V[M], OV[N];
    for (int j=0; j<M; j++)
      V[j] = B[j][i];

    matXvec<N,M,OP>(A,V,OV);

    for (int j=0; j<N; j++)
      O[j][i] = OV[j];
  }
}


CCS_MAIN(int, char**) {
  
  const int N=8, M=8, K=4;

  T A[N][M], B[M][K], V[M], OV[N], OM[N][K];

  float HI = 6546364.0;
  float LO = -6546364.0;
  srand(time(NULL));
  
  for (int j=0; j<M; j++) {
    V[j] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO))); 
    for (int i=0; i<N; i++) 
      A[i][j] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));

    for (int i=0; i<K; i++)
      B[j][i] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
  }

  const OP_TYPE OP = SEPARATE;
  // const OP_TYPE OP = FMA;
  // const OP_TYPE OP = DOT;

  matXvec<N,M,OP>(A,V,OV);

  matXmat<N,M,K,OP>(A,B,OM);

  std::cout << "MATRIX x VECTOR\n" << std::endl;
  for (int i=0; i<N; i++)
    std::cout << OV[i].to_float() << std::endl;
  
  std::cout << "\nMATRIX x MATRIX\n" << std::endl;
  for (int i=0; i<N; i++){
    for (int j=0; j<K; j++) {
      if (OM[i][j].sign)
        std::cout << OM[i][j].to_float() << "   ";
      else
        std::cout << OM[i][j].to_float() << "    ";
    }
    std::cout << std::endl;
  }
}