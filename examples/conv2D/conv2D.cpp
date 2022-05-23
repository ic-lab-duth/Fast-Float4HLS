// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build an 2D Convolution
#include <iostream>
#include <random>
#include "fast_float.h"
#include "ac_channel.h"

#include "mc_scverify.h"

typedef ffp16b T;
typedef ac_channel<T> chT;

#pragma hls_design
template<int IF_HEIGHT,int IF_WIDTH, int KERNEL>
void CCS_BLOCK(conv2D) (chT &inp, T filter[KERNEL][KERNEL], chT &out) {
  
  static T lb[KERNEL-1][IF_WIDTH];
  static T window[KERNEL][KERNEL];

  #pragma hls_pipeline_init_interval 1
  for (int i=0; i<IF_HEIGHT; i++) {
    for (int j=0; j<IF_WIDTH; j++) {

      T cur_pxl = inp.read();
      
      // Window and Line buffers update
      #pragma hls_unroll
      for (int m=0; m<KERNEL; m++) {
        T tmp = (m < KERNEL-1) ? lb[m][j] : cur_pxl;
        #pragma hls_unroll
        for (int n=0; n<KERNEL; n++) {
          window[m][n] = (n<KERNEL-1) ? window[m][n+1] : tmp;
        }
        if (m>0)
          lb[m-1][j] = tmp;
      }
      
      // Output Pixel Calculation
      if (i>=KERNEL-1 && j>=KERNEL-1) {
        T o_pxl = 0.0;
        #pragma hls_unroll
        for (int m=0; m<KERNEL; m++) {
          T tmp_pxl;
          tmp_pxl.dotProd<KERNEL>(window[m], filter[m]);
          o_pxl += tmp_pxl;
        }
        out.write(o_pxl);
      }
    }
  }
}

float rand_num(float HI, float LO) {
  return  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
}

CCS_MAIN(int, char**) {
  
  const int H = 14;
  const int W = 14;
  const int K = 5;
  
  chT img, out;
  T flt[K][K];

  float HI;
  float LO;

  srand(time(NULL));

  for (int i=0; i<H; i++)
    for (int j=0; j<W; j++) {
      HI = rand();
      LO = -rand();
      img.write(rand_num(HI,LO));
    }

  for (int i=0; i<K; i++)
    for (int j=0; j<K; j++) {
      HI = rand();
      LO = -rand();
      flt[i][j] = rand_num(HI,LO);
    }

  conv2D<H,W,K>(img,flt,out);


  for (int i=0; i<H-K+1; i++) {
    for (int j=0; j<W-K+1; j++) 
      std::cout << out.read().to_float() << " ";
    
    std::cout << std::endl;
  }

  CCS_RETURN(0);
}
