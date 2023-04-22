//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include <fmt/core.h>
#include <gtest/gtest.h>

#include <cmath>

#define EPS 1e-8

extern "C" {
#include "star/matrix4.h"
}

TEST(Matrix4, VectorMultiplication) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat x[4] = {1, 2, 3, 4};
  sfloat y[4];
  star_VecMul44(y, A, x);
  EXPECT_EQ(y[0], 90);
  EXPECT_EQ(y[1], 100);
  EXPECT_EQ(y[2], 110);
  EXPECT_EQ(y[3], 120);

  sfloat At[16];
  star_Transpose44(At, A);
  star_VecMul44(y, At, x);
  EXPECT_EQ(y[0], 30);
  EXPECT_EQ(y[1], 70);
  EXPECT_EQ(y[2], 110);
  EXPECT_EQ(y[3], 150);

  star_TransposedVecMul44(y, A, x);
  EXPECT_EQ(y[0], 30);
  EXPECT_EQ(y[1], 70);
  EXPECT_EQ(y[2], 110);
  EXPECT_EQ(y[3], 150);
}

TEST(Matrix4, MatrixMultiplication) {
  // Normal multiplication
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat B[16] = {0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8,
                  0.9, -1.0, 1.1, -1.2, 1.3, -1.4, 1.5, -1.6};
  sfloat C[16];
  sfloat C_expected[16] = {-3.4, -3.6,  -3.8,  -4.0,  -6.6,  -6.8,  -7.0,  -7.2,
                           -9.8, -10.0, -10.2, -10.4, -13.0, -13.2, -13.4, -13.6};
  star_MatMul44(C, A, B);

  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], C_expected[i], EPS);
  }

  // Check transposes
  sfloat At[16];
  sfloat Bt[16];
  sfloat Ct[16];
  star_Transpose44(At, A);
  star_Transpose44(Bt, B);
  star_Transpose44(Ct, C_expected);
  star_MatMul44(C, Bt, At);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], Ct[i], EPS);
  }

  // C = At*B
  star_MatMul44(C_expected, At, B);
  star_TransposedMatMul44(C, A, B);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], C_expected[i], EPS);
  }

  // C = A*Bt
  star_MatMul44(C_expected, A, Bt);
  star_MatMulTransposed44(C, A, B);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], C_expected[i], EPS);
  }
}

TEST(Matrix4, SetZero) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  star_SetZero44(A);
  for (int i = 0; i < 16; i++) {
    EXPECT_EQ(A[i], 0);
  }
}

TEST(Matrix4, SetConst) {
  sfloat A[16];
  star_SetConst44(A, 3.14);
  for (int i = 0; i < 16; i++) {
    EXPECT_EQ(A[i], 3.14);
  }
}

TEST(Matrix4, SetIdentity) {
  sfloat A[16];
  star_SetIdentity44(A);
  for (int i = 0; i < 16; i++) {
    if (i % 5 == 0) {
      EXPECT_EQ(A[i], 1);
    } else {
      EXPECT_EQ(A[i], 0);
    }
  }
}

TEST(Matrix4, SetDiagonal) {
  sfloat A[16];
  sfloat d[4] = {1, 2, 3, 4};
  star_SetDiagonal44(A, d);
  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(A[i + 4 * i], d[i], EPS);
  }
}

TEST(Matrix4, Copy) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat B[16];
  star_Copy44(B, A);
  for (int i = 0; i < 16; i++) {
    EXPECT_EQ(B[i], A[i]);
  }
}

TEST(Matrix4, Add) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat B[16] = {0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8,
                  0.9, -1.0, 1.1, -1.2, 1.3, -1.4, 1.5, -1.6};
  sfloat C[16];
  star_Add44(C, A, B);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] + B[i], EPS);
  }
}

TEST(Matrix4, Sub) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat B[16] = {0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8,
                  0.9, -1.0, 1.1, -1.2, 1.3, -1.4, 1.5, -1.6};
  sfloat C[16];
  star_Sub44(C, A, B);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] - B[i], EPS);
  }
}

TEST(Matrix4, Mul) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat B[16] = {0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8,
                  0.9, -1.0, 1.1, -1.2, 1.3, -1.4, 1.5, -1.6};
  sfloat C[16];
  star_Mul44(C, A, B);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] * B[i], EPS);
  }
}

TEST(Matrix4, Div) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat B[16] = {0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8,
                  0.9, -1.0, 1.1, -1.2, 1.3, -1.4, 1.5, -1.6};
  sfloat C[16];
  star_Div44(C, A, B);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] / B[i], EPS);
  }
}

TEST(Matrix4, AddConst) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat C[16];
  sfloat alpha = 0.5;
  star_AddConst44(C, A, alpha);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] + alpha, EPS);
  }
}

TEST(Matrix4, SubConst) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat C[16];
  sfloat alpha = 0.5;
  star_SubConst44(C, A, alpha);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] - alpha, EPS);
  }
}

TEST(Matrix4, MulConst) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat C[16];
  sfloat alpha = 0.5;
  star_MulConst44(C, A, alpha);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] * alpha, EPS);
  }
}

TEST(Matrix4, DivConst) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat C[16];
  sfloat alpha = 0.5;
  star_DivConst44(C, A, alpha);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(C[i], A[i] / alpha, EPS);
  }
}

TEST(Matrix4, AddConst_Aliased) {
  sfloat A[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  sfloat A0[16];
  star_Copy44(A0, A);
  sfloat alpha = 0.5;
  star_AddConst44(A, A, alpha);
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR(A[i], A0[i] + alpha, EPS);
  }
}
