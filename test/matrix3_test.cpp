//
// Created by Brian Jackson on 4/2/23.
// Copyright (c) 2023. All rights reserved.
//

#include <cmath>
#include <gtest/gtest.h>

extern "C" {
#include "star/matrix3.h"
}

TEST(Matrix3, SetZero) {
  sfloat mat[9];
  star_SetZero33(mat);
  EXPECT_EQ(mat[0], 0);
  EXPECT_EQ(mat[1], 0);
  EXPECT_EQ(mat[2], 0);
  EXPECT_EQ(mat[3], 0);
  EXPECT_EQ(mat[4], 0);
  EXPECT_EQ(mat[5], 0);
  EXPECT_EQ(mat[6], 0);
  EXPECT_EQ(mat[7], 0);
  EXPECT_EQ(mat[8], 0);
}

TEST(Matrix3, SetConst) {
  sfloat mat[9];
  const sfloat value = 1.2;
  star_SetConst33(mat, value);
  EXPECT_EQ(mat[0], value);
  EXPECT_EQ(mat[1], value);
  EXPECT_EQ(mat[2], value);
  EXPECT_EQ(mat[3], value);
  EXPECT_EQ(mat[4], value);
  EXPECT_EQ(mat[5], value);
  EXPECT_EQ(mat[6], value);
  EXPECT_EQ(mat[7], value);
  EXPECT_EQ(mat[8], value);
}

TEST(Matrix3, SetIdentity) {
  sfloat mat[9];
  const sfloat value = 1.2;
  star_SetIdentity33(mat, value);
  EXPECT_EQ(mat[0], value);
  EXPECT_EQ(mat[1], 0);
  EXPECT_EQ(mat[2], 0);
  EXPECT_EQ(mat[3], 0);
  EXPECT_EQ(mat[4], value);
  EXPECT_EQ(mat[5], 0);
  EXPECT_EQ(mat[6], 0);
  EXPECT_EQ(mat[7], 0);
  EXPECT_EQ(mat[8], value);
}

TEST(Matrix3, SetDiagonal) {
  sfloat mat[9];
  star_SetZero33(mat);
  const sfloat diag[3] = {1.2, 2.3, 3.4};
  star_SetDiagonal33(mat, diag);
  EXPECT_EQ(mat[0], diag[0]);
  EXPECT_EQ(mat[1], 0);
  EXPECT_EQ(mat[2], 0);
  EXPECT_EQ(mat[3], 0);
  EXPECT_EQ(mat[4], diag[1]);
  EXPECT_EQ(mat[5], 0);
  EXPECT_EQ(mat[6], 0);
  EXPECT_EQ(mat[7], 0);
  EXPECT_EQ(mat[8], diag[2]);
}

TEST(Matrix3, MulAB) {
  sfloat C[9];
  const sfloat A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const sfloat B[9] = {9, 8, 7, 6, 5, 4, 3, 2, 1};
  star_MatMul33(C, A, B);
  sfloat C_expected[9] = {90, 114, 138, 54, 69, 84, 18, 24, 30};
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }
}

TEST(Matrix3, MulAtB) {
  sfloat C[9];
  const sfloat At[9] = {1, 4, 7, 2, 5, 8, 3, 6, 9};
  const sfloat B[9] = {9, 8, 7, 6, 5, 4, 3, 2, 1};
  star_TransposedMatMul33(C, At, B);
  sfloat C_expected[9] = {90, 114, 138, 54, 69, 84, 18, 24, 30};
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }
}

TEST(Matrix3, MulABt) {
  sfloat C[9];
  const sfloat A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const sfloat Bt[9] = {9, 6, 3, 8, 5, 2, 7, 4, 1};
  star_MatMulTransposed33(C, A, Bt);
  sfloat C_expected[9] = {90, 114, 138, 54, 69, 84, 18, 24, 30};
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }
}

TEST(Matrix3, MulAx) {
  sfloat y[3];
  const sfloat A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const sfloat x[3] = {10, 11, 12};
  star_VecMul33(y, A, x);
  sfloat y_expected[3] = {138, 171, 204};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(y[i], y_expected[i]);
  }
}

TEST(Matrix3, MulAtx) {
  sfloat y[3];
  const sfloat At[9] = {1, 4, 7, 2, 5, 8, 3, 6, 9};
  const sfloat x[3] = {10, 11, 12};
  star_TransposedVecMul33(y, At, x);
  sfloat y_expected[3] = {138, 171, 204};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(y[i], y_expected[i]);
  }
}

TEST(Matrix3, MatMul_UpperTriangular) {
  sfloat A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  sfloat U[9] = {1, 0, 0, 2, 3, 0, 4, 5, 6};
  sfloat C[9];
  star_UpperMatMul33(C, U, A);
  sfloat C_expected[9] = {17, 21, 18, 38, 45, 36, 59, 69, 54};
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }

  // set lower triangular portion of U to nonzero values to make sure it is the same
  U[1] = 1;
  U[2] = 2;
  U[5] = 3;
  star_UpperMatMul33(C, U, A);
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }

  // test with a vector
  sfloat x[3] = {1, 2, 3};
  sfloat y[3];
  star_UpperVecMul33(y, U, x);
  sfloat y_expected[3] = {17, 21, 18};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(y[i], y_expected[i]);
  }
}

TEST(Matrix3, MatMul_LowerTriangular) {
  sfloat A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  sfloat L[9] = {1, 2, 3, 0, 4, 5, 0, 0, 6};
  sfloat C[9];
  star_LowerMatMul33(C, L, A);
  sfloat C_expected[9] = {1, 10, 31, 4, 28, 73, 7, 46, 115};
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }

  // set upper triangular portion to nonzero and make sure it is ignored
  L[3] = 10;
  L[6] = 11;
  L[7] = 12;
  star_LowerMatMul33(C, L, A);
  for (int i = 0; i < 9; i++) {
    EXPECT_EQ(C[i], C_expected[i]);
  }

  // test with vector
  sfloat x[3] = {1, 2, 3};
  sfloat y[3];
  star_LowerVecMul33(y, L, x);
  sfloat y_expected[3] = {1, 10, 31};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(y[i], y_expected[i]);
  }
}

TEST(Matrix3, UpperTriSolve) {
  sfloat U[9] = {1, 0, 0, 2, 3, 0, 4, 5, 6};
  sfloat b[3] = {14, 13, 12};
  sfloat x[3];
  star_UpperTriSolve33(x, U, b);
  sfloat x_expected[3] = {4, 1, 2};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(x[i], x_expected[i]);
  }

  // set lower triangular portion of U to nonzero values to make sure it is the same
  U[1] = 1;
  U[2] = 2;
  U[5] = 3;
  star_UpperTriSolve33(x, U, b);
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(x[i], x_expected[i]);
  }

  // Check U*x = b
  sfloat b2[3];
  star_UpperVecMul33(b2, U, x);
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(b2[i], b[i]);
  }
}

TEST(Matrix3, LowerTriSolve) {
  sfloat L[9] = {1, 2, 3, 0, 4, 5, 0, 0, 6};
  sfloat b[3] = {2, -4, 14};
  sfloat x[3];
  star_LowerTriSolve33(x, L, b);
  sfloat x_expected[3] = {2, -2, 3};
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(x[i], x_expected[i]);
  }

  // set upper triangular portion to nonzero and make sure it is ignored
  L[3] = 10;
  L[6] = 11;
  L[7] = 12;
  star_LowerTriSolve33(x, L, b);
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(x[i], x_expected[i]);
  }

  // Check L*x = b
  sfloat b2[3];
  star_LowerVecMul33(b2, L, x);
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(b2[i], b[i]);
  }
}

TEST(Matrix3, Determinant) {
  sfloat A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  sfloat det = star_Det33(A);
  EXPECT_EQ(det, 0);
  A[0] = 2;
  det = star_Det33(A);
  EXPECT_EQ(det, -3);

  // Check determinant for a rotation matrix
  sfloat theta = 0.5;
  sfloat R[9] = {cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1};
  det = star_Det33(R);
  EXPECT_EQ(det, 1);
}