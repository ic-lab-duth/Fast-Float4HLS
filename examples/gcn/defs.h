// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network
#ifndef _DEFS_H
#define _DEFS_H

#include "fast_float.h"
#include "ac_int.h"
#include "matrix.h"

typedef ffp16b btype;


/*****Sizes for Citeseer *****/
// size of the problem (matrix dimensions)
static const int N = 3327;
static const int I_F = 3703;
static const int O_F1 = 21;
static const int O_F2 = 6;

// number of non-zero elements of the sparse matrix
static const int nZ = 12431;

// number of non-zero elements of the feature matrix
static const int feat_nZ = 105165;

/*****Sizes for Cora *****/
/*// size of the problem (matrix dimensions)
static const int N = 2708;
static const int I_F = 1433;
static const int O_F1 = 64;
static const int O_F2 = 7;

// number of non-zero elements of the sparse matrix
static const int nZ = 13264;

// number of non-zero elements of the feature matrix
static const int feat_nZ = 49216;*/

/*****Sizes for Pubmed *****/
/*// size of the problem (matrix dimensions)
static const int N = 19717;
static const int I_F = 500;
static const int O_F1 = 18;
static const int O_F2 = 3;

// number of non-zero elements of the sparse matrix
static const int nZ = 108365;

// number of non-zero elements of the feature matrix
static const int feat_nZ = 988031;*/

#endif
