// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network 
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>

#include "defs.h"
#include "helper.h"
#include "gcn.h"
#include "mc_scverify.h"


CCS_MAIN(int argc, char** argv) {
  
  gcn< btype, N, I_F, O_F1, O_F2, nZ > GCN;

  btype* h_o = new btype[N*O_F2];
  btype* h = new btype[N*I_F];

  btype* w1 = new btype[I_F*O_F1];
  btype* w2 = new btype[O_F1*O_F2];
  
  ac_int<32, false>* a_row = new ac_int<32, false>[N+1];
  ac_int<32, false>* a_col = new ac_int<32, false>[nZ];
  btype* a_val = new btype[nZ];
  
  int* A_row = new int[N+1];
  int* A_col = new int[nZ];
  float* A_val = new float[nZ];

  int* H_row = new int[N+1];
  int* H_col = new int[feat_nZ];
  float* H_val = new float[feat_nZ];

  float* W1 = new float[I_F*O_F1];
  float* W2 = new float[O_F1*O_F2];
  
  float* H = new float[N*I_F];


  // read input matrices from txt files
  read_csr<float, N, nZ>(A_row, A_col, A_val, "./matrices/citeseer_adj.txt");
  read_csr<float, N, feat_nZ>(H_row, H_col, H_val, "./matrices/citeseer_feat.txt");
  read_data<float, I_F, O_F1>(W1, "./matrices/citeseer_weights.txt");
  read_data<float, O_F1, O_F2>(W2, "./matrices/citeseer_weights2.txt");

  to_array<float, N, I_F, feat_nZ>(H_row, H_col, H_val, H);

  for (int i=0; i < N+1; i++) {
    a_row[i] = A_row[i];
  }

  for (int i=0; i < nZ; i++) {
    a_col[i] = A_col[i];
    a_val[i] = A_val[i];
  }

  for (int i=0; i < N*I_F; i++) {
    h[i] = H[i];
  }

  for (int i=0; i < I_F*O_F1; i++) {
    w1[i] = W1[i];
  }

  for (int i=0; i < O_F1*O_F2; i++) {
    w2[i] = W2[i];
  }

  //initialize output
  for (int i=0; i < N; i++) {
    for (int j=0; j < O_F2; j++) {
      h_o[i*O_F2 + j] = 0.0;
    }
  }

  GCN.run_wrap(a_row, a_col, a_val, h, w1, w2, h_o);

  std::cout << "DONE!!!" << std::endl;

  CCS_RETURN(0);
}
