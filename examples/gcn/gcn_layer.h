// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network
#ifndef _GCN_LAYER_H
#define _GCN_LAYER_H

#include "defs.h"
#include "mc_scverify.h"

#pragma hls_design
template< typename T, int N, int M, int K , int L, int nZ, int TILE_HEIGHT, int TILE_WIDTH>
class gcn_layer{
private:
    
    // ReLU activation function
    void relu (T a[L], T b[L]) {
        COL: for (int j = 0; j < L; j++) {
            b[j] = (a[j] > 0) ? a[j] : 0;
        }
    }
    
    // In the future will be replaced by softmax
    // Currently implement only hard decision of the most popular class 
    void argmax (T in[L], T out[L]) {
        T maxVal = -65535.0;
        T curVal;
        MAX_ITER: for (int j=0; j < L; j++) {
            curVal = in[j];
            if (curVal > maxVal) {
                maxVal = curVal;
            }
        }
                
        COL: for (int j=0; j < L; j++) {
            out[j] = (in[j]==maxVal) ? (T)1.0 : (T)0.0;
        }
    }

public:
    gcn_layer(){}

    #pragma hls_design
    void CCS_BLOCK(run) (T a_val, T h_row[TILE_HEIGHT], T w[TILE_HEIGHT][TILE_WIDTH], T h_o[L], short int p, short int q, bool change_row, bool last_layer, bool apply_activation) {

      static T ah_prod[TILE_HEIGHT];
      static T out_buf[L];
      
      // reset buffer when row changed
      if(change_row) {
        for (int i=0; i < L; i++) {
          out_buf[i] = 0.0;
        }
      }

      for (int i=0; i < TILE_HEIGHT; i++) {
        ah_prod[i] = a_val * h_row[i];
      }

      for (int i=0; i < TILE_WIDTH; i++) {
        for (int j=0; j < TILE_HEIGHT; j++) {
          if (q*TILE_WIDTH + i < L)
            out_buf[q*TILE_WIDTH + i] += w[j][i] * ah_prod[j];
        }
      }

      // when computation of an output row completed apply activation function
      // we assume the use of ReLU activation for the hidden layers and hard decision for the output layer
      if (apply_activation) {
        if (!last_layer) {
          relu(out_buf, h_o);
        } else {
          argmax(out_buf, h_o);
        }
      }
      
    }

};

#endif
