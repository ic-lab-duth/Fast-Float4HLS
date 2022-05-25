// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network
#ifndef _GCN_H
#define _GCN_H

#include "defs.h"
#include "mc_scverify.h"

#pragma hls_design
template< typename T, int N, int I_F, int O_F1, int O_F2, int nZ >
class gcn{
private:
    
    // ReLU activation function
    void relu (T a[O_F1], T b[O_F1]) {
        COL: for (int j = 0; j < O_F1; j++) {
            b[j] = (a[j] > 0) ? a[j] : 0;
        }
    }
    
    // In the future will be replaced by softmax
    // Currently implement only hard decision of the most popular class 
    void argmax (T in[O_F2], T out[O_F2]) {
        T maxVal = -65535.0;
        T curVal;
        MAX_ITER: for (int j=0; j < O_F2; j++) {
            curVal = in[j];
            if (curVal > maxVal) {
                maxVal = curVal;
            }
        }
                
        COL: for (int j=0; j < O_F2; j++) {
            out[j] = (in[j]==maxVal) ? (T)1.0 : (T)0.0;
        }
    }

    void run (Array< ac_int<32,false> > a_row,
              Array< ac_int<32,false> > a_col,
              Array<T> a_val,
              Matrix<T> h_i,
              Matrix<T> w1,
              Matrix<T> w2,
              Matrix<T> h_o) {

      #ifndef __SYNTHESIS__
      T* inter_buf = new T[N*O_F1];
      Matrix<T> inter(N, O_F1, inter_buf);
      #else
      static T inter_buf[N*O_F1];
      Matrix<T> inter(N, O_F1, inter_buf);
      #endif

      // first GCN layer
      ac_int<32,false> cur_row = a_row[0][0];
      ac_int<32,false> next_row;

      static T out_buf1[O_F1];

      A_ROW: for (int i=0; i < N; i++) {
        if ( i < 100 ) {

        // initialize output buffer
        for (int j=0; j < O_F1; j++) {
          out_buf1[j] = 0.0;
        }

        next_row = a_row[0][i+1];
        ROW_NZ: for (int j=cur_row; j < next_row; j++) {
          INP_FEAT: for (int p=0; p < I_F; p++) {
            OUT_FEAT: for (int q=0; q < O_F1; q++) {
              out_buf1[q] += (a_val[0][j] * h_i[a_col[0][j]][p]) * w1[p][q];
            }
          }
        }
      

        relu(out_buf1, out_buf1);
        WR_OUT:for (int j=0; j < O_F1; j++) {
          // write the output
          inter[i][j] = out_buf1[j];
        }

        cur_row = next_row;
        }
      }

      // second GCN layer

      cur_row = a_row[0][0];
      
      static T out_buf2[O_F2];

      A_ROW_2: for (int i=0; i < N; i++) {
        
        if ( i < 100 ) {
        // initialize output buffer
        for (int j=0; j < O_F1; j++) {
          out_buf2[j] = 0.0;
        }

        next_row = a_row[0][i+1];
        ROW_NZ_2: for (int j=cur_row; j < next_row; j++) {
          INP_FEAT_2: for (int p=0; p < O_F1; p++) {
            OUT_FEAT_2: for (int q=0; q < O_F2; q++) {
              out_buf2[q] += (a_val[0][j] * inter[a_col[0][j]][p]) * w2[p][q];
            }
          }
        }
      

        argmax(out_buf2, out_buf2);
        WR_OUT_2:for (int j=0; j < O_F1; j++) {
          // write the output
          h_o[i][j] = out_buf2[j];
        }

        cur_row = next_row;
        }
      }

    }

public:
    gcn(){}

    #pragma hls_design interface
    void CCS_BLOCK(run_wrap) (ac_int<32, false> a_row[N+1],
                              ac_int<32, false> a_col[nZ],
                              T a_val[nZ],
                              T h_i[N*I_F],
                              T w1[I_F*O_F1],
                              T w2[I_F*O_F2],
                              T h_o[N*O_F2]){

      Array< ac_int<32, false> > a_row_(N+1, a_row);
      Array< ac_int<32, false> > a_col_(nZ, a_col);
      Array< T > a_val_(nZ, a_val);
      Matrix< T > h_i_(N, I_F, h_i);
      Matrix< T > h_o_(N, O_F2, h_o);
      Matrix< T > w1_(I_F, O_F1, w1);
      Matrix< T > w2_(O_F1, O_F2, w2);
      
      run(a_row_, a_col_, a_val_, h_i_, w1_, w2_, h_o_);
    }
    

};

#endif
