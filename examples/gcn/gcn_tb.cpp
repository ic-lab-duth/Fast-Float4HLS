// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network
#include <iostream>

#include "defs.h"
#include "helper.h"
#include "gcn_layer.h"
#include "mc_scverify.h"


CCS_MAIN(int argc, char** argv) {
  
  // number of vertical and horizontal tiles for first layer
  const int ver1 = ((I_F1 % K)==0) ? (I_F1 / K) : (I_F1 / K)+1;
  const int hor1 = ((O_F1 % M)==0) ? (O_F1 / M) : (O_F1 / M)+1;
  
  // number of vertical and horizontal tiles for second layer
  const int ver2 = ((O_F1 % K)==0) ? (O_F1 / K) : (O_F1 / K)+1;
  const int hor2 = ((O_F2 % M)==0) ? (O_F2 / M) : (O_F2 / M)+1;

  gcn_layer< btype, N,N,I_F1,O_F1,nZ, K, M> layer1;
  gcn_layer< btype, N,N,O_F1,O_F2,nZ, K, M> layer2;

  btype h_o[N][O_F2];
  btype inter_feats[N][O_F1];

  std::vector<int> a_row;
  std::vector<int> a_col;
  std::vector<float> a_val;
  std::vector< std::vector<float> > h;
  std::vector< std::vector<float> > W1;
  std::vector< std::vector<float> > W2;

  // read input matrices from txt files for citeseer
  read_adj<float, N, nZ>(a_row, a_col, a_val, "./matrices/citeseer_adj.txt");
  read_data<float, N, I_F1>(h, "./matrices/citeseer_feat.txt");
  read_data<float, I_F1, O_F1>(W1, "./matrices/citeseer_weights.txt");
  read_data<float, O_F1, O_F2>(W2, "./matrices/citeseer_weights2.txt");

/*  // read input matrices from txt files for cora
  read_adj<float, N, nZ>(a_row, a_col, a_val, "./matrices/cora_adj.txt");
  read_data<float, N, I_F1>(h, "./matrices/cora_feat.txt");
  read_data<float, I_F1, O_F1>(W1, "./matrices/cora_weights.txt");
  read_data<float, O_F1, O_F2>(W2, "./matrices/cora_weights2.txt");*/

/*  // read input matrices from txt files for pubmed
  read_adj<float, N, nZ>(a_row, a_col, a_val, "./matrices/pubmed_adj.txt");
  read_data<float, N, I_F1>(h, "./matrices/pubmed_feat.txt");
  read_data<float, I_F1, O_F1>(W1, "./matrices/pubmed_weights.txt");
  read_data<float, O_F1, O_F2>(W2, "./matrices/pubmed_weights2.txt");*/

  
  btype w_tile[K][M];


  // first gcn layer
  std::cout << "Executing first GCN layer..." << std::endl;
  bool change_row = false;
  bool last_layer = false;
  bool apply_activation = false;
  int lambda;
  btype a_lambda;
  btype out_row1[O_F1];
  btype in_row1[K];

  for (int i=0; i < N; i++) {
    change_row = true;  // flag that we change row in adjacency matrix

    for (int j=a_row[i]; j < a_row[i+1]; j++) {
      lambda = a_col[j];
      a_lambda = a_val[j];


      for (int p=0; p < ver1; p++) {
        for (int r=0; r < K; r++) {
          if (p*K+r < I_F1)
            in_row1[r] = h[lambda][p*K + r];
        }
        
        for (int q=0; q < hor1; q++) {
          
          // read the weight tile used in this step
          for (int r=0; r < K; r++) {
            if ((p*K+r) < I_F1) {
              for (int s=0; s < M; s++) {
                if ((q*M+s) < O_F1) {
                  w_tile[r][s] = W1.at(p*K+r)[q*M+s];
                }
              }
            }
          }
          
          // when the computation of an output row completed, apply activation function
          if(j==a_row[i+1]-1)
            apply_activation = true;

          layer1.run(a_lambda, in_row1, w_tile, out_row1, p, q, change_row, last_layer, apply_activation);
          // reset flags
          change_row = false;
          apply_activation = false;
        }
      }

    }
    
    WR_OUT_INTER: for (int j=0; j < O_F1; j++) {
      // write the output
      inter_feats[i][j] = out_row1[j];
    }
  }


  // second gcn layer
  std::cout << "Executing second GCN layer..." << std::endl;
  last_layer = true;
  btype out_row2[O_F2];
  btype in_row2[K];

  for (int i=0; i < N; i++) {
    change_row = true;  // flag that we change row in adjacency matrix
    
    for (int j=a_row[i]; j < a_row[i+1]; j++) {
      lambda = a_col[j];
      a_lambda = a_val[j];


      for (int p=0; p < ver2; p++) {
        for (int r=0; r < K; r++) {
          if (p*K+r < O_F1)
            in_row2[r] = h[lambda][p*K + r];
        }
        
        for (int q=0; q < hor2; q++) {
          
          // read the weight tile used in this step
          for (int r=0; r < K; r++) {
            if ((p*K+r) < O_F1) {
              for (int s=0; s < M; s++) {
                if ((q*M+s) < O_F2) {
                  w_tile[r][s] = W2.at(p*K+r)[q*M+s];
                }
              }
            }
          }
          
          // when the computation of an output row completed, apply activation function
          if(j==a_row[i+1]-1)
            apply_activation = true;

          layer2.run(a_lambda, in_row2, w_tile, out_row2, p, q, change_row, last_layer, apply_activation);
          // reset flags
          change_row = false;
          apply_activation = false;
        }
      }

    }
    
    WR_OUT_FINAL: for (int j=0; j < O_F2; j++) {
      // write the output
      h_o[i][j] = out_row2[j];
    }

  }

  std::cout << "DONE!!!" << std::endl;

  CCS_RETURN(0);
}
