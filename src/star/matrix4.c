//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "matrix4.h"

#define IDX(i, j) ((i) + (j)*4)

void star_SetZero44(sfloat* mat) {
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
  mat[12] = 0;
  mat[13] = 0;
  mat[14] = 0;
  mat[15] = 0;
}

void star_SetConst44(sfloat* mat, sfloat value) {
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
  mat[12] = value;
  mat[13] = value;
  mat[14] = value;
  mat[15] = value;
}

void star_SetIdentity44(sfloat* mat) {
  mat[0] = 1;
  mat[1] = 0;
  mat[2] = 0;
  mat[3] = 0;
  mat[4] = 0;
  mat[5] = 1;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = 0;
  mat[9] = 0;
  mat[10] = 1;
  mat[11] = 0;
  mat[12] = 0;
  mat[13] = 0;
  mat[14] = 0;
  mat[15] = 1;
}

void star_SetDiagonal44(sfloat* mat, const sfloat* diag) {
  mat[0] = diag[0];
  mat[5] = diag[1];
  mat[10] = diag[2];
  mat[15] = diag[3];
}

void star_Copy44(sfloat* dst, const sfloat* src) {
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
  dst[12] = src[12];
  dst[13] = src[13];
  dst[14] = src[14];
  dst[15] = src[15];
}

void star_Transpose44(sfloat* dst, const sfloat* src) {
  dst[IDX(0, 0)] = src[IDX(0, 0)];
  dst[IDX(0, 1)] = src[IDX(1, 0)];
  dst[IDX(0, 2)] = src[IDX(2, 0)];
  dst[IDX(0, 3)] = src[IDX(3, 0)];
  dst[IDX(1, 0)] = src[IDX(0, 1)];
  dst[IDX(1, 1)] = src[IDX(1, 1)];
  dst[IDX(1, 2)] = src[IDX(2, 1)];
  dst[IDX(1, 3)] = src[IDX(3, 1)];
  dst[IDX(2, 0)] = src[IDX(0, 2)];
  dst[IDX(2, 1)] = src[IDX(1, 2)];
  dst[IDX(2, 2)] = src[IDX(2, 2)];
  dst[IDX(2, 3)] = src[IDX(3, 2)];
  dst[IDX(3, 0)] = src[IDX(0, 3)];
  dst[IDX(3, 1)] = src[IDX(1, 3)];
  dst[IDX(3, 2)] = src[IDX(2, 3)];
  dst[IDX(3, 3)] = src[IDX(3, 3)];
}

void star_TransposeInPlace44(sfloat* mat) {
  sfloat tmp;
  tmp = mat[IDX(1, 0)];
  mat[IDX(1, 0)] = mat[IDX(0, 1)];
  mat[IDX(0, 1)] = tmp;
  tmp = mat[IDX(2, 0)];
  mat[IDX(2, 0)] = mat[IDX(0, 2)];
  mat[IDX(0, 2)] = tmp;
  tmp = mat[IDX(3, 0)];
  mat[IDX(3, 0)] = mat[IDX(0, 3)];
  mat[IDX(0, 3)] = tmp;
  tmp = mat[IDX(2, 1)];
  mat[IDX(2, 1)] = mat[IDX(1, 2)];
  mat[IDX(1, 2)] = tmp;
  tmp = mat[IDX(3, 1)];
  mat[IDX(3, 1)] = mat[IDX(1, 3)];
  mat[IDX(1, 3)] = tmp;
  tmp = mat[IDX(3, 2)];
  mat[IDX(3, 2)] = mat[IDX(2, 3)];
  mat[IDX(2, 3)] = tmp;
}

void star_MatMul44(sfloat* C, const sfloat* A, const sfloat* B) {
  // C = A * B, where A, B, and C are 4x4 matrices stored column-major
  C[0] = A[0] * B[0] + A[4] * B[1] + A[8] * B[2] + A[12] * B[3];
  C[1] = A[1] * B[0] + A[5] * B[1] + A[9] * B[2] + A[13] * B[3];
  C[2] = A[2] * B[0] + A[6] * B[1] + A[10] * B[2] + A[14] * B[3];
  C[3] = A[3] * B[0] + A[7] * B[1] + A[11] * B[2] + A[15] * B[3];
  C[4] = A[0] * B[4] + A[4] * B[5] + A[8] * B[6] + A[12] * B[7];
  C[5] = A[1] * B[4] + A[5] * B[5] + A[9] * B[6] + A[13] * B[7];
  C[6] = A[2] * B[4] + A[6] * B[5] + A[10] * B[6] + A[14] * B[7];
  C[7] = A[3] * B[4] + A[7] * B[5] + A[11] * B[6] + A[15] * B[7];
  C[8] = A[0] * B[8] + A[4] * B[9] + A[8] * B[10] + A[12] * B[11];
  C[9] = A[1] * B[8] + A[5] * B[9] + A[9] * B[10] + A[13] * B[11];
  C[10] = A[2] * B[8] + A[6] * B[9] + A[10] * B[10] + A[14] * B[11];
  C[11] = A[3] * B[8] + A[7] * B[9] + A[11] * B[10] + A[15] * B[11];
  C[12] = A[0] * B[12] + A[4] * B[13] + A[8] * B[14] + A[12] * B[15];
  C[13] = A[1] * B[12] + A[5] * B[13] + A[9] * B[14] + A[13] * B[15];
  C[14] = A[2] * B[12] + A[6] * B[13] + A[10] * B[14] + A[14] * B[15];
  C[15] = A[3] * B[12] + A[7] * B[13] + A[11] * B[14] + A[15] * B[15];
}

void star_VecMul44(sfloat* y, const sfloat* A, const sfloat* x) {
  // Extract out x so that x and y can be aliased
  sfloat x0 = x[0];
  sfloat x1 = x[1];
  sfloat x2 = x[2];
  sfloat x3 = x[3];

  y[0] = A[0] * x0 + A[4] * x1 + A[8] * x2 + A[12] * x3;
  y[1] = A[1] * x0 + A[5] * x1 + A[9] * x2 + A[13] * x3;
  y[2] = A[2] * x0 + A[6] * x1 + A[10] * x2 + A[14] * x3;
  y[3] = A[3] * x0 + A[7] * x1 + A[11] * x2 + A[15] * x3;
}
void star_TransposedMatMul44(sfloat* C, const sfloat* At, const sfloat* B) {
  // C = A^T * B, where A, B, and C are 4x4 matrices stored column-major
  C[0] = At[0] * B[0] + At[1] * B[1] + At[2] * B[2] + At[3] * B[3];
  C[1] = At[4] * B[0] + At[5] * B[1] + At[6] * B[2] + At[7] * B[3];
  C[2] = At[8] * B[0] + At[9] * B[1] + At[10] * B[2] + At[11] * B[3];
  C[3] = At[12] * B[0] + At[13] * B[1] + At[14] * B[2] + At[15] * B[3];

  C[4] = At[0] * B[4] + At[1] * B[5] + At[2] * B[6] + At[3] * B[7];
  C[5] = At[4] * B[4] + At[5] * B[5] + At[6] * B[6] + At[7] * B[7];
  C[6] = At[8] * B[4] + At[9] * B[5] + At[10] * B[6] + At[11] * B[7];
  C[7] = At[12] * B[4] + At[13] * B[5] + At[14] * B[6] + At[15] * B[7];

  C[8] = At[0] * B[8] + At[1] * B[9] + At[2] * B[10] + At[3] * B[11];
  C[9] = At[4] * B[8] + At[5] * B[9] + At[6] * B[10] + At[7] * B[11];
  C[10] = At[8] * B[8] + At[9] * B[9] + At[10] * B[10] + At[11] * B[11];
  C[11] = At[12] * B[8] + At[13] * B[9] + At[14] * B[10] + At[15] * B[11];

  C[12] = At[0] * B[12] + At[1] * B[13] + At[2] * B[14] + At[3] * B[15];
  C[13] = At[4] * B[12] + At[5] * B[13] + At[6] * B[14] + At[7] * B[15];
  C[14] = At[8] * B[12] + At[9] * B[13] + At[10] * B[14] + At[11] * B[15];
  C[15] = At[12] * B[12] + At[13] * B[13] + At[14] * B[14] + At[15] * B[15];
}

void star_MatMulTransposed44(sfloat* C, const sfloat* A, const sfloat* Bt) {
  // C = A * B^T, where A, B, and C are 4x4 matrices stored column-major
  C[0] = A[0] * Bt[0] + A[4] * Bt[4] + A[8] * Bt[8] + A[12] * Bt[12];
  C[1] = A[1] * Bt[0] + A[5] * Bt[4] + A[9] * Bt[8] + A[13] * Bt[12];
  C[2] = A[2] * Bt[0] + A[6] * Bt[4] + A[10] * Bt[8] + A[14] * Bt[12];
  C[3] = A[3] * Bt[0] + A[7] * Bt[4] + A[11] * Bt[8] + A[15] * Bt[12];

  C[4] = A[0] * Bt[1] + A[4] * Bt[5] + A[8] * Bt[9] + A[12] * Bt[13];
  C[5] = A[1] * Bt[1] + A[5] * Bt[5] + A[9] * Bt[9] + A[13] * Bt[13];
  C[6] = A[2] * Bt[1] + A[6] * Bt[5] + A[10] * Bt[9] + A[14] * Bt[13];
  C[7] = A[3] * Bt[1] + A[7] * Bt[5] + A[11] * Bt[9] + A[15] * Bt[13];

  C[8] = A[0] * Bt[2] + A[4] * Bt[6] + A[8] * Bt[10] + A[12] * Bt[14];
  C[9] = A[1] * Bt[2] + A[5] * Bt[6] + A[9] * Bt[10] + A[13] * Bt[14];
  C[10] = A[2] * Bt[2] + A[6] * Bt[6] + A[10] * Bt[10] + A[14] * Bt[14];
  C[11] = A[3] * Bt[2] + A[7] * Bt[6] + A[11] * Bt[10] + A[15] * Bt[14];

  C[12] = A[0] * Bt[3] + A[4] * Bt[7] + A[8] * Bt[11] + A[12] * Bt[15];
  C[13] = A[1] * Bt[3] + A[5] * Bt[7] + A[9] * Bt[11] + A[13] * Bt[15];
  C[14] = A[2] * Bt[3] + A[6] * Bt[7] + A[10] * Bt[11] + A[14] * Bt[15];
  C[15] = A[3] * Bt[3] + A[7] * Bt[7] + A[11] * Bt[11] + A[15] * Bt[15];
}

void star_TransposedVecMul44(sfloat* y, const sfloat* At, const sfloat* x) {
  // Extract out x so that x and y can be aliased
  sfloat x0 = x[0];
  sfloat x1 = x[1];
  sfloat x2 = x[2];
  sfloat x3 = x[3];

  // y = A^T * x, A stored column-major, not transposed
  y[0] = At[0] * x0 + At[1] * x1 + At[2] * x2 + At[3] * x3;
  y[1] = At[4] * x0 + At[5] * x1 + At[6] * x2 + At[7] * x3;
  y[2] = At[8] * x0 + At[9] * x1 + At[10] * x2 + At[11] * x3;
  y[3] = At[12] * x0 + At[13] * x1 + At[14] * x2 + At[15] * x3;
}

void star_Add44(sfloat* C, const sfloat* A, const sfloat* B) {
  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];
  C[4] = A[4] + B[4];
  C[5] = A[5] + B[5];
  C[6] = A[6] + B[6];
  C[7] = A[7] + B[7];
  C[8] = A[8] + B[8];
  C[9] = A[9] + B[9];
  C[10] = A[10] + B[10];
  C[11] = A[11] + B[11];
  C[12] = A[12] + B[12];
  C[13] = A[13] + B[13];
  C[14] = A[14] + B[14];
  C[15] = A[15] + B[15];
}

void star_Sub44(sfloat* C, const sfloat* A, const sfloat* B) {
  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];
  C[4] = A[4] - B[4];
  C[5] = A[5] - B[5];
  C[6] = A[6] - B[6];
  C[7] = A[7] - B[7];
  C[8] = A[8] - B[8];
  C[9] = A[9] - B[9];
  C[10] = A[10] - B[10];
  C[11] = A[11] - B[11];
  C[12] = A[12] - B[12];
  C[13] = A[13] - B[13];
  C[14] = A[14] - B[14];
  C[15] = A[15] - B[15];
}

void star_Mul44(sfloat* C, const sfloat* A, const sfloat* B) {
  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];
  C[4] = A[4] * B[4];
  C[5] = A[5] * B[5];
  C[6] = A[6] * B[6];
  C[7] = A[7] * B[7];
  C[8] = A[8] * B[8];
  C[9] = A[9] * B[9];
  C[10] = A[10] * B[10];
  C[11] = A[11] * B[11];
  C[12] = A[12] * B[12];
  C[13] = A[13] * B[13];
  C[14] = A[14] * B[14];
  C[15] = A[15] * B[15];
}

void star_Div44(sfloat* C, const sfloat* A, const sfloat* B) {
  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];
  C[4] = A[4] / B[4];
  C[5] = A[5] / B[5];
  C[6] = A[6] / B[6];
  C[7] = A[7] / B[7];
  C[8] = A[8] / B[8];
  C[9] = A[9] / B[9];
  C[10] = A[10] / B[10];
  C[11] = A[11] / B[11];
  C[12] = A[12] / B[12];
  C[13] = A[13] / B[13];
  C[14] = A[14] / B[14];
  C[15] = A[15] / B[15];
}

void star_AddConst44(sfloat* C, const sfloat* A, sfloat b) {
  C[0] = A[0] + b;
  C[1] = A[1] + b;
  C[2] = A[2] + b;
  C[3] = A[3] + b;
  C[4] = A[4] + b;
  C[5] = A[5] + b;
  C[6] = A[6] + b;
  C[7] = A[7] + b;
  C[8] = A[8] + b;
  C[9] = A[9] + b;
  C[10] = A[10] + b;
  C[11] = A[11] + b;
  C[12] = A[12] + b;
  C[13] = A[13] + b;
  C[14] = A[14] + b;
  C[15] = A[15] + b;
}

void star_SubConst44(sfloat* C, const sfloat* A, sfloat b) {
  C[0] = A[0] - b;
  C[1] = A[1] - b;
  C[2] = A[2] - b;
  C[3] = A[3] - b;
  C[4] = A[4] - b;
  C[5] = A[5] - b;
  C[6] = A[6] - b;
  C[7] = A[7] - b;
  C[8] = A[8] - b;
  C[9] = A[9] - b;
  C[10] = A[10] - b;
  C[11] = A[11] - b;
  C[12] = A[12] - b;
  C[13] = A[13] - b;
  C[14] = A[14] - b;
  C[15] = A[15] - b;
}

void star_MulConst44(sfloat* C, const sfloat* A, sfloat b) {
  C[0] = A[0] * b;
  C[1] = A[1] * b;
  C[2] = A[2] * b;
  C[3] = A[3] * b;
  C[4] = A[4] * b;
  C[5] = A[5] * b;
  C[6] = A[6] * b;
  C[7] = A[7] * b;
  C[8] = A[8] * b;
  C[9] = A[9] * b;
  C[10] = A[10] * b;
  C[11] = A[11] * b;
  C[12] = A[12] * b;
  C[13] = A[13] * b;
  C[14] = A[14] * b;
  C[15] = A[15] * b;
}

void star_DivConst44(sfloat* C, const sfloat* A, sfloat b) {
  C[0] = A[0] / b;
  C[1] = A[1] / b;
  C[2] = A[2] / b;
  C[3] = A[3] / b;
  C[4] = A[4] / b;
  C[5] = A[5] / b;
  C[6] = A[6] / b;
  C[7] = A[7] / b;
  C[8] = A[8] / b;
  C[9] = A[9] / b;
  C[10] = A[10] / b;
  C[11] = A[11] / b;
  C[12] = A[12] / b;
  C[13] = A[13] / b;
  C[14] = A[14] / b;
  C[15] = A[15] / b;
}
