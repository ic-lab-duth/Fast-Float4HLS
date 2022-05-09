// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Example using Fast-Float4HLS to build a Graph Convolutional Network
#ifndef _DEFS_H
#define _DEFS_H

#include "fast_float.h"
#include "mc_scverify.h"

typedef float dtype;
typedef ffp32 btype;


/*****Sizes for Citeseer *****/
// size of the problem (matrix dimensions)
static const int N = 3327;
static const int I_F1 = 3703;
static const int O_F1 = 21;
static const int O_F2 = 6;

// number of non-zero elements of the sparse matrix
static const int nZ = 12431;

/*****Sizes for Cora *****/
/*// size of the problem (matrix dimensions)
static const int N = 2708;
static const int I_F1 = 1433;
static const int O_F1 = 64;
static const int O_F2 = 7;

// number of non-zero elements of the sparse matrix
static const int nZ = 13264;*/

/*****Sizes for Pubmed *****/
/*// size of the problem (matrix dimensions)
static const int N = 19717;
static const int I_F1 = 500;
static const int O_F1 = 18;
static const int O_F2 = 3;

// number of non-zero elements of the sparse matrix
static const int nZ = 108365;*/


// tile size
static const int K = 16;  // rows
static const int M = 16;  // columns


#endif
