//
// Created by Brian Jackson on 3/9/23.
// Copyright (c) 2023. All rights reserved.
//

#include "matrix3.h"

void star_SetZero33(sfloat mat[9]) {
  mat[0] = 0;
  mat[1] = 0;
  mat[2] = 0;
  mat[3] = 0;
  mat[4] = 0;
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = 0;
}

void star_SetConst33(sfloat mat[9], sfloat value) {
  mat[0] = value;
  mat[1] = value;
  mat[2] = value;
  mat[3] = value;
  mat[4] = value;
  mat[5] = value;
  mat[6] = value;
  mat[7] = value;
  mat[8] = value;
}

void star_SetIdentity33(sfloat mat[9], sfloat value) {
  mat[0] = value;
  mat[1] = 0;
  mat[2] = 0;
  mat[3] = 0;
  mat[4] = value;
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = value;
}

void star_SetDiagonal33(sfloat mat[9], const sfloat diag[3]) {
  mat[0] = diag[0];
  mat[1] = 0;
  mat[2] = 0;
  mat[3] = 0;
  mat[4] = diag[1];
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = diag[2];
}

void star_MatMul33(sfloat* C, const sfloat* A, const sfloat* B) {
  C[0] = A[0] * B[0] + A[3] * B[1] + A[6] * B[2];
  C[1] = A[1] * B[0] + A[4] * B[1] + A[7] * B[2];
  C[2] = A[2] * B[0] + A[5] * B[1] + A[8] * B[2];
  C[3] = A[0] * B[3] + A[3] * B[4] + A[6] * B[5];
  C[4] = A[1] * B[3] + A[4] * B[4] + A[7] * B[5];
  C[5] = A[2] * B[3] + A[5] * B[4] + A[8] * B[5];
  C[6] = A[0] * B[6] + A[3] * B[7] + A[6] * B[8];
  C[7] = A[1] * B[6] + A[4] * B[7] + A[7] * B[8];
  C[8] = A[2] * B[6] + A[5] * B[7] + A[8] * B[8];
}

void star_MatTransposedMul33(sfloat* C, const sfloat* At, const sfloat* B) {
  C[0] = At[0] * B[0] + At[1] * B[1] + At[2] * B[2];
  C[1] = At[3] * B[0] + At[4] * B[1] + At[5] * B[2];
  C[2] = At[6] * B[0] + At[7] * B[1] + At[8] * B[2];
  C[3] = At[0] * B[3] + At[1] * B[4] + At[2] * B[5];
  C[4] = At[3] * B[3] + At[4] * B[4] + At[5] * B[5];
  C[5] = At[6] * B[3] + At[7] * B[4] + At[8] * B[5];
  C[6] = At[0] * B[6] + At[1] * B[7] + At[2] * B[8];
  C[7] = At[3] * B[6] + At[4] * B[7] + At[5] * B[8];
  C[8] = At[6] * B[6] + At[7] * B[7] + At[8] * B[8];
}

void star_MatMulTransposed33(sfloat* C, const sfloat* A, const sfloat* Bt) {
  C[0] = A[0] * Bt[0] + A[3] * Bt[3] + A[6] * Bt[6];
  C[1] = A[1] * Bt[0] + A[4] * Bt[3] + A[7] * Bt[6];
  C[2] = A[2] * Bt[0] + A[5] * Bt[3] + A[8] * Bt[6];
  C[3] = A[0] * Bt[1] + A[3] * Bt[4] + A[6] * Bt[7];
  C[4] = A[1] * Bt[1] + A[4] * Bt[4] + A[7] * Bt[7];
  C[5] = A[2] * Bt[1] + A[5] * Bt[4] + A[8] * Bt[7];
  C[6] = A[0] * Bt[2] + A[3] * Bt[5] + A[6] * Bt[8];
  C[7] = A[1] * Bt[2] + A[4] * Bt[5] + A[7] * Bt[8];
  C[8] = A[2] * Bt[2] + A[5] * Bt[5] + A[8] * Bt[8];
}

void star_VecMul33(sfloat* C, const sfloat* A, const sfloat* x) {
  C[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  C[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  C[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

void star_VecTransposedMul33(sfloat* C, const sfloat* At, const sfloat* x) {
  C[0] = At[0] * x[0] + At[1] * x[1] + At[2] * x[2];
  C[1] = At[3] * x[0] + At[4] * x[1] + At[5] * x[2];
  C[2] = At[6] * x[0] + At[7] * x[1] + At[8] * x[2];
}

void star_UpperMatMul33(sfloat* C, const sfloat* U, const sfloat* A) {
  C[0] = U[0] * A[0] + U[3] * A[1] + U[6] * A[2];
  C[1] = U[4] * A[1] + U[7] * A[2];
  C[2] = U[8] * A[2];
  C[3] = U[0] * A[3] + U[3] * A[4] + U[6] * A[5];
  C[4] = U[4] * A[4] + U[7] * A[5];
  C[5] = U[8] * A[5];
  C[6] = U[0] * A[6] + U[3] * A[7] + U[6] * A[8];
  C[7] = U[4] * A[7] + U[7] * A[8];
  C[8] = U[8] * A[8];
}

void star_UpperVecMul33(sfloat* C, const sfloat* U, const sfloat* x) {
  C[0] = U[0] * x[0] + U[3] * x[1] + U[6] * x[2];
  C[1] = U[4] * x[1] + U[7] * x[2];
  C[2] = U[8] * x[2];
}

void star_LowerMatMul33(sfloat* C, const sfloat* L, const sfloat* A) {
  C[0] = L[0] * A[0];
  C[1] = L[1] * A[0] + L[4] * A[1];
  C[2] = L[2] * A[0] + L[5] * A[1] + L[8] * A[2];
  C[3] = L[0] * A[3];
  C[4] = L[1] * A[3] + L[4] * A[4];
  C[5] = L[2] * A[3] + L[5] * A[4] + L[8] * A[5];
  C[6] = L[0] * A[6];
  C[7] = L[1] * A[6] + L[4] * A[7];
  C[8] = L[2] * A[6] + L[5] * A[7] + L[8] * A[8];
}

void star_LowerVecMul33(sfloat* C, const sfloat* L, const sfloat* x) {
  C[0] = L[0] * x[0];
  C[1] = L[1] * x[0] + L[4] * x[1];
  C[2] = L[2] * x[0] + L[5] * x[1] + L[8] * x[2];
}

void star_UpperTriSolve33(sfloat* x, const sfloat* U, const sfloat* b) {
  x[2] = b[2] / U[8];
  x[1] = (b[1] - U[7] * x[2]) / U[4];
  x[0] = (b[0] - U[3] * x[1] - U[6] * x[2]) / U[0];
}

void star_LowerTriSolve33(sfloat* x, const sfloat* L, const sfloat* b) {
  x[0] = b[0] / L[0];
  x[1] = (b[1] - L[1] * x[0]) / L[4];
  x[2] = (b[2] - L[2] * x[0] - L[5] * x[1]) / L[8];
}

sfloat star_Det33(const sfloat* mat) {
  return mat[0] * mat[4] * mat[8] + mat[3] * mat[7] * mat[2] + mat[6] * mat[1] * mat[5] -
         mat[6] * mat[4] * mat[2] - mat[0] * mat[7] * mat[5] - mat[3] * mat[1] * mat[8];
}
