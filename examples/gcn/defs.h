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

// size of the problem (matrix dimensions)
static const int N = 3327;
static const int I_F1 = 3703;
static const int O_F1 = 21;
static const int O_F2 = 6;

// tile sizes
static const int K = 16;  // rows
static const int M = 16;  // columns

// number of non-zero elements of the sparse matrix
static const int nZ = 12431;


#endif
