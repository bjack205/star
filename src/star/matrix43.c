//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "matrix43.h"

#define IDX(i, j) ((i) + 4 * (j))

void star_SetZero43(sfloat mat[12]) {
  mat[0] = 0;
  mat[1] = 0;
  mat[2] = 0;
  mat[3] = 0;
  mat[4] = 0;
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = 0;
  mat[9] = 0;
  mat[10] = 0;
  mat[11] = 0;
}

void star_SetConst43(sfloat mat[12], sfloat value) {
  mat[0] = value;
  mat[1] = value;
  mat[2] = value;
  mat[3] = value;
  mat[4] = value;
  mat[5] = value;
  mat[6] = value;
  mat[7] = value;
  mat[8] = value;
  mat[9] = value;
  mat[10] = value;
  mat[11] = value;
}

void star_Copy43(sfloat dst[12], const sfloat src[12]) {
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];
  dst[3] = src[3];
  dst[4] = src[4];
  dst[5] = src[5];
  dst[6] = src[6];
  dst[7] = src[7];
  dst[8] = src[8];
  dst[9] = src[9];
  dst[10] = src[10];
  dst[11] = src[11];
}

//void star_Transpose43(sfloat dst[12], const sfloat src[12]) {
//  dst[IDX(0, 0)] = src[IDX(0, 0)];
//  dst[IDX(0, 1)] = src[IDX(1, 0)];
//  dst[IDX(0, 2)] = src[IDX(2, 0)];
//  dst[IDX(0, 3)] = src[IDX(3, 0)];
//  dst[IDX(1, 0)] = src[IDX(0, 1)];
//  dst[IDX(1, 1)] = src[IDX(1, 1)];
//  dst[IDX(1, 2)] = src[IDX(2, 1)];
//  dst[IDX(1, 3)] = src[IDX(3, 1)];
//  dst[IDX(2, 0)] = src[IDX(0, 2)];
//  dst[IDX(2, 1)] = src[IDX(1, 2)];
//  dst[IDX(2, 2)] = src[IDX(2, 2)];
//  dst[IDX(2, 3)] = src[IDX(3, 2)];
//}

/*---------------------------------*/
/* Multiplication                  */
/*---------------------------------*/

void star_MatMul433(sfloat C43[12], const sfloat A43[12], const sfloat B33[9]) {
  // Multiply a 4x3 matrix by a 3x3 matrix
  C43[0] = A43[0] * B33[0] + A43[4] * B33[1] + A43[8] * B33[2];
  C43[1] = A43[1] * B33[0] + A43[5] * B33[1] + A43[9] * B33[2];
  C43[2] = A43[2] * B33[0] + A43[6] * B33[1] + A43[10] * B33[2];
  C43[3] = A43[3] * B33[0] + A43[7] * B33[1] + A43[11] * B33[2];

  C43[4] = A43[0] * B33[3] + A43[4] * B33[4] + A43[8] * B33[5];
  C43[5] = A43[1] * B33[3] + A43[5] * B33[4] + A43[9] * B33[5];
  C43[6] = A43[2] * B33[3] + A43[6] * B33[4] + A43[10] * B33[5];
  C43[7] = A43[3] * B33[3] + A43[7] * B33[4] + A43[11] * B33[5];

  C43[8] = A43[0] * B33[6] + A43[4] * B33[7] + A43[8] * B33[8];
  C43[9] = A43[1] * B33[6] + A43[5] * B33[7] + A43[9] * B33[8];
  C43[10] = A43[2] * B33[6] + A43[6] * B33[7] + A43[10] * B33[8];
  C43[11] = A43[3] * B33[6] + A43[7] * B33[7] + A43[11] * B33[8];
}

void star_MatMulTransposed433(sfloat C43[12], const sfloat A43[12], const sfloat B33t[9]) {
  C43[0] = A43[0] * B33t[0] + A43[4] * B33t[3] + A43[8] * B33t[6];
  C43[1] = A43[1] * B33t[0] + A43[5] * B33t[3] + A43[9] * B33t[6];
  C43[2] = A43[2] * B33t[0] + A43[6] * B33t[3] + A43[10] * B33t[6];
  C43[3] = A43[3] * B33t[0] + A43[7] * B33t[3] + A43[11] * B33t[6];

  C43[4] = A43[0] * B33t[1] + A43[4] * B33t[4] + A43[8] * B33t[7];
  C43[5] = A43[1] * B33t[1] + A43[5] * B33t[4] + A43[9] * B33t[7];
  C43[6] = A43[2] * B33t[1] + A43[6] * B33t[4] + A43[10] * B33t[7];
  C43[7] = A43[3] * B33t[1] + A43[7] * B33t[4] + A43[11] * B33t[7];

  C43[8] = A43[0] * B33t[2] + A43[4] * B33t[5] + A43[8] * B33t[8];
  C43[9] = A43[1] * B33t[2] + A43[5] * B33t[5] + A43[9] * B33t[8];
  C43[10] = A43[2] * B33t[2] + A43[6] * B33t[5] + A43[10] * B33t[8];
  C43[11] = A43[3] * B33t[2] + A43[7] * B33t[5] + A43[11] * B33t[8];
}

void star_MatMul443(sfloat C43[12], const sfloat A44[16], const sfloat B43[12]) {
  // Multiply a 4x4 matrix by a 4x3 matrix
  C43[0] = A44[0] * B43[0] + A44[4] * B43[1] + A44[8] * B43[2] + A44[12] * B43[3];
  C43[1] = A44[1] * B43[0] + A44[5] * B43[1] + A44[9] * B43[2] + A44[13] * B43[3];
  C43[2] = A44[2] * B43[0] + A44[6] * B43[1] + A44[10] * B43[2] + A44[14] * B43[3];
  C43[3] = A44[3] * B43[0] + A44[7] * B43[1] + A44[11] * B43[2] + A44[15] * B43[3];

  C43[4] = A44[0] * B43[4] + A44[4] * B43[5] + A44[8] * B43[6] + A44[12] * B43[7];
  C43[5] = A44[1] * B43[4] + A44[5] * B43[5] + A44[9] * B43[6] + A44[13] * B43[7];
  C43[6] = A44[2] * B43[4] + A44[6] * B43[5] + A44[10] * B43[6] + A44[14] * B43[7];
  C43[7] = A44[3] * B43[4] + A44[7] * B43[5] + A44[11] * B43[6] + A44[15] * B43[7];

  C43[8] = A44[0] * B43[8] + A44[4] * B43[9] + A44[8] * B43[10] + A44[12] * B43[11];
  C43[9] = A44[1] * B43[8] + A44[5] * B43[9] + A44[9] * B43[10] + A44[13] * B43[11];
  C43[10] = A44[2] * B43[8] + A44[6] * B43[9] + A44[10] * B43[10] + A44[14] * B43[11];
  C43[11] = A44[3] * B43[8] + A44[7] * B43[9] + A44[11] * B43[10] + A44[15] * B43[11];
}

void star_TransposedMatMul443(sfloat C43[12], const sfloat A44t[16], const sfloat B43[12]) {
  C43[0] = A44t[0] * B43[0] + A44t[1] * B43[1] + A44t[2] * B43[2] + A44t[3] * B43[3];
  C43[1] = A44t[4] * B43[0] + A44t[5] * B43[1] + A44t[6] * B43[2] + A44t[7] * B43[3];
  C43[2] = A44t[8] * B43[0] + A44t[9] * B43[1] + A44t[10] * B43[2] + A44t[11] * B43[3];
  C43[3] = A44t[12] * B43[0] + A44t[13] * B43[1] + A44t[14] * B43[2] + A44t[15] * B43[3];

  C43[4] = A44t[0] * B43[4] + A44t[1] * B43[5] + A44t[2] * B43[6] + A44t[3] * B43[7];
  C43[5] = A44t[4] * B43[4] + A44t[5] * B43[5] + A44t[6] * B43[6] + A44t[7] * B43[7];
  C43[6] = A44t[8] * B43[4] + A44t[9] * B43[5] + A44t[10] * B43[6] + A44t[11] * B43[7];
  C43[7] = A44t[12] * B43[4] + A44t[13] * B43[5] + A44t[14] * B43[6] + A44t[15] * B43[7];

  C43[8] = A44t[0] * B43[8] + A44t[1] * B43[9] + A44t[2] * B43[10] + A44t[3] * B43[11];
  C43[9] = A44t[4] * B43[8] + A44t[5] * B43[9] + A44t[6] * B43[10] + A44t[7] * B43[11];
  C43[10] = A44t[8] * B43[8] + A44t[9] * B43[9] + A44t[10] * B43[10] + A44t[11] * B43[11];
  C43[11] = A44t[12] * B43[8] + A44t[13] * B43[9] + A44t[14] * B43[10] + A44t[15] * B43[11];
}

void star_MatMul344(sfloat C34[12], const sfloat A34[12], const sfloat B44[16]) {
  (void)C34;
  (void)A34;
  (void)B44;
}

void star_MatMulTransposed344(sfloat C34[12], const sfloat A34[12], const sfloat B44t[16]) {
  (void)C34;
  (void)A34;
  (void)B44t;
}

void star_VecMul43(sfloat y[4], const sfloat A[12], const sfloat x[3]) {
  // Multiply a 4x3 matrix by a 3-vector
  y[0] = A[0] * x[0] + A[4] * x[1] + A[8] * x[2];
  y[1] = A[1] * x[0] + A[5] * x[1] + A[9] * x[2];
  y[2] = A[2] * x[0] + A[6] * x[1] + A[10] * x[2];
  y[3] = A[3] * x[0] + A[7] * x[1] + A[11] * x[2];
}

void star_TransposedVecMul43(sfloat y[3], const sfloat At[12], const sfloat x[4]) {
  // Multiply a 3x4 matrix by a 4-vector
  y[0] = At[0] * x[0] + At[1] * x[1] + At[2] * x[2] + At[3] * x[3];
  y[1] = At[4] * x[0] + At[5] * x[1] + At[6] * x[2] + At[7] * x[3];
  y[2] = At[8] * x[0] + At[9] * x[1] + At[10] * x[2] + At[11] * x[3];
}
