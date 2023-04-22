//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include <fmt/core.h>
#include <gtest/gtest.h>

#include "star/Mat3.hpp"
#include "star/Transpose.hpp"
#include "star/Vec3.hpp"
#include "star/matrix_multiplication.hpp"

#define EPS 1e-8

using namespace star;

TEST(Matrix3, Zero) {
  Mat3 A = Mat3::Zero();
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], 0, EPS);
  }
}

TEST(Matrix3, InitializerList) {
  Mat3 A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], i + 1, EPS);
  }
}

TEST(Matrix3, ByRows) {
  // clang-format off
  Mat3 A = Mat3::ByRows(
      1, 4, 7,
      2, 5, 8,
      3, 6, 9
  );
  // clang-format on
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], i + 1, EPS);
  }
}

TEST(Matrix3, VectorMultiplication) {
  Mat3 A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec3 v = {1, 2, 3};
  Vec3 v2 = A * v;
  EXPECT_NEAR(v2[0], 30, EPS);
  EXPECT_NEAR(v2[1], 36, EPS);
  EXPECT_NEAR(v2[2], 42, EPS);
}

TEST(Matrix3, TransposedVectorMultiplication) {
  Mat3 A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Vec3 x = {1, 2, 3};
  Vec3 y0 = Multiply(Transpose(A), x);
  EXPECT_NEAR(y0[0], 14, EPS);
  EXPECT_NEAR(y0[1], 32, EPS);
  EXPECT_NEAR(y0[2], 50, EPS);

  Vec3 y1 = Transpose(A) * x;
  EXPECT_NEAR(y1[0], 14, EPS);
  EXPECT_NEAR(y1[1], 32, EPS);
  EXPECT_NEAR(y1[2], 50, EPS);
}