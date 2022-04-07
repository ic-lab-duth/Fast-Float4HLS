#include "fast_float.h"
#include "mc_scverify.h"


template<int N, typename T = ffp32>
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
}

CCS_MAIN(int, char**) {
  
  typedef ffp32 DATA_T;

  const int NSIZE = 10;

  DATA_T inA[NSIZE];

  float HI = 65415456364.0;
  float LO = -654364.0;

  srand(time(NULL));

  for (int i=0; i<NSIZE; i++) {
    inA[i] =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
    std::cout << inA[i].to_float() << std::endl;
  }
  std::cout << std::endl;
  
  bubbleSort<NSIZE>(inA);

  for (int i=0; i<NSIZE; i++) {
    std::cout << inA[i].to_float() << std::endl;
  }

  CCS_RETURN(0);
}