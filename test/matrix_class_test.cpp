//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include <fmt/core.h>
#include <gtest/gtest.h>

#include "star/Mat3.hpp"
#include "star/Mat4.hpp"
#include "star/Mat43.hpp"
#include "star/Transpose.hpp"
#include "star/Vec3.hpp"
#include "star/Vec4.hpp"
#include "star/matrix_multiplication.hpp"

#define EPS 1e-8

using namespace star;

TEST(Matrix3, Zero) {
  Mat3 A = Mat3::Zero();
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], 0, EPS);
  }
  EXPECT_EQ(A.Rows(), 3);
  EXPECT_EQ(A.Cols(), 3);
  EXPECT_EQ(A.Size(), 9);
  EXPECT_EQ(Mat3::kRows, 3);
  EXPECT_EQ(Mat3::kCols, 3);
  EXPECT_EQ(Mat3::kSize, 9);
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

TEST(Matrix4, Zero) {
  Mat4 A = Mat4::Zero();
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], 0, EPS);
  }
  EXPECT_EQ(A.Rows(), 4);
  EXPECT_EQ(A.Cols(), 4);
  EXPECT_EQ(A.Size(), 16);
  EXPECT_EQ(Mat4::kRows, 4);
  EXPECT_EQ(Mat4::kCols, 4);
  EXPECT_EQ(Mat4::kSize, 16);
}

TEST(Matrix4, InitializerList) {
  Mat4 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], i + 1, EPS);
  }

  // Check that it's stored column-major
  EXPECT_NEAR(A(1, 0), 2, EPS);
  EXPECT_NEAR(A(0, 1), 5, EPS);
}

TEST(Matrix4, ByRows) {
  // clang-format off
  Mat4 A = Mat4::ByRows(
    1, 5, 9, 13,
    2, 6, 10, 14,
    3, 7, 11, 15,
    4, 8, 12, 16
  );
  // clang-format on

  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], i + 1, EPS);
  }

  // Check that it's stored column-major
  EXPECT_NEAR(A(1, 0), 2, EPS);
  EXPECT_NEAR(A(0, 1), 5, EPS);
}

TEST(Matrix4, VectorMultiplication) {
  Mat4 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  Vec4 x = {1, 2, 3, 4};
  Vec4 y = A * x;
  EXPECT_NEAR(y[0], 90, EPS);
  EXPECT_NEAR(y[1], 100, EPS);
  EXPECT_NEAR(y[2], 110, EPS);
  EXPECT_NEAR(y[3], 120, EPS);
}

TEST(Matrix4, TransposedVectorMultiplication) {
  Mat4 At = Mat4::ByRows(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
  Vec4 x = {1, 2, 3, 4};
  Vec4 y = Multiply(Transpose(At), x);
  EXPECT_NEAR(y[0], 90, EPS);
  EXPECT_NEAR(y[1], 100, EPS);
  EXPECT_NEAR(y[2], 110, EPS);
  EXPECT_NEAR(y[3], 120, EPS);

  Mat4 A = At.Transpose();
  Vec4 y2 = A * x;
  EXPECT_NEAR(y2[0], 90, EPS);
  EXPECT_NEAR(y2[1], 100, EPS);
  EXPECT_NEAR(y2[2], 110, EPS);
  EXPECT_NEAR(y2[3], 120, EPS);
}

TEST(Matrix4, MatrixMultiplication) {
  Mat4 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  Mat4 B = Mat4::ByRows(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
  Mat4 C_expected = {276, 304, 332, 360, 304, 336, 368, 400,
                     332, 368, 404, 440, 360, 400, 440, 480};
  Mat4 C = A * B;
  for (int i = 0; i < C.Size(); ++i) {
    EXPECT_NEAR(C[i], C_expected[i], EPS);
  }
}

TEST(Matrix4, TransposedMatrixMultiplication) {
  Mat4 A = {-8, 3, -5, -1, -10, -6, 9, -7, -1, -2, -7, -4, 10, -9, 9, 2};
  Mat4 B = {-4, -2, -10, -9, 0, 5, -1, -2, 1, -10, -8, 8, -5, -10, -7, 8};
  Mat4 At = A.Transpose();
  Mat4 Bt = B.Transpose();

  // Check A'B
  Mat4 C1 = At * B;
  Mat4 C2 = Transpose(A) * B;
  for (int i = 0; i < C1.Size(); ++i) {
    EXPECT_NEAR(C1[i], C2[i], EPS);
  }

  // Check A*B'
  C1 = A * Bt;
  C2 = A * Transpose(B);
  for (int i = 0; i < C1.Size(); ++i) {
    EXPECT_NEAR(C1[i], C2[i], EPS);
  }

  // Check C' = B'*A'
  Mat4 Ct = (A * B).Transpose();
  C1 = Bt * At;
  for (int i = 0; i < C1.Size(); ++i) {
    EXPECT_NEAR(C1[i], Ct[i], EPS);
  }
}

TEST(Matrix43, Zero) {
  Mat43 A = Mat43::Zero();
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], 0, EPS);
  }
  EXPECT_EQ(A.Rows(), 4);
  EXPECT_EQ(A.Cols(), 3);
  EXPECT_EQ(A.Size(), 12);
  EXPECT_EQ(Mat43::kRows, 4);
  EXPECT_EQ(Mat43::kCols, 3);
  EXPECT_EQ(Mat43::kSize, 12);
}

TEST(Matrix43, InitializerList) {
  Mat43 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], i + 1, EPS);
  }

  // Check that it's stored column-major
  EXPECT_NEAR(A(1, 0), 2, EPS);
  EXPECT_NEAR(A(0, 1), 5, EPS);
}

TEST(Matrix43, Indexing) {
  Mat43 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  for (int i = 0; i < A.Rows(); ++i) {
    for (int j = 0; j < A.Cols(); ++j) {
      EXPECT_NEAR(A(i, j), i + 4 * j + 1, EPS);
      sfloat val = A[{i, j}];
      EXPECT_NEAR(val, i + 4 * j + 1, EPS);
      EXPECT_NEAR(A[i + 4 * j], i + 4 * j + 1, EPS);
    }
  }
}

TEST(Matrix43, ByRows) {
  // clang-format off
  Mat43 A = Mat43::ByRows(
      1, 5, 9,
      2, 6, 10,
      3, 7, 11,
      4, 8, 12
  );
  // clang-format on
  for (int i = 0; i < A.Size(); ++i) {
    EXPECT_NEAR(A[i], i + 1, EPS);
  }

  // Check that it's stored column-major
  EXPECT_NEAR(A(0, 1), 5, EPS);
  EXPECT_NEAR(A(1, 0), 2, EPS);
}

TEST(Matrix43, VectorMultiplication) {
  Mat43 A = Mat43::ByRows(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
  Vec3 x = {1, 2, 3};
  Vec4 y = Multiply(A, x);
  EXPECT_NEAR(y[0], 14, EPS);
  EXPECT_NEAR(y[1], 32, EPS);
  EXPECT_NEAR(y[2], 50, EPS);
  EXPECT_NEAR(y[3], 68, EPS);

  Vec4 y2 = A * x;
  EXPECT_NEAR(y2[0], 14, EPS);
  EXPECT_NEAR(y2[1], 32, EPS);
  EXPECT_NEAR(y2[2], 50, EPS);
  EXPECT_NEAR(y2[3], 68, EPS);
}

TEST(Matrix43, TransposedVectorMultiplication) {
  Mat43 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Vec4 x = {1, 2, 3, 4};
  Vec3 y = Multiply(Transpose(A), x);
  EXPECT_NEAR(y[0], 30, EPS);
  EXPECT_NEAR(y[1], 70, EPS);
  EXPECT_NEAR(y[2], 110, EPS);

  Vec3 y2 = Transpose(A) * x;
  EXPECT_NEAR(y2[0], 30, EPS);
  EXPECT_NEAR(y2[1], 70, EPS);
  EXPECT_NEAR(y2[2], 110, EPS);

  Vec3 y3 = Vec3::Zero();
  MultiplyInPlace(y3, Transpose(A), x);
}

TEST(Matrix43, MatrixMultiplication) {
  Mat43 A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Mat4 B = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  Mat3 C = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  Mat43 BA_expected = {90, 100, 110, 120, 202, 228, 254, 280, 314, 356, 398, 440};
  Mat43 BA = B * A;
  for (int i = 0; i < BA.Size(); ++i) {
    EXPECT_NEAR(BA[i], BA_expected[i], EPS);
  }

  Mat43 AC_expected = {38, 44, 50, 56, 83, 98, 113, 128, 128, 152, 176, 200};
  Mat43 AC = A * C;
  for (int i = 0; i < AC.Size(); ++i) {
    EXPECT_NEAR(AC[i], AC_expected[i], EPS);
  }

  // Transpose the square matrices
  Mat43 BtA_expected = {30, 70, 110, 150, 70, 174, 278, 382, 110, 278, 446, 614};
  Mat43 BtA = Transpose(B) * A;
  for (int i = 0; i < BtA.Size(); ++i) {
    EXPECT_NEAR(BtA[i], BtA_expected[i], EPS);
  }

  Mat43 ACt_expected = {84, 96, 108, 120, 99, 114, 129, 144, 114, 132, 150, 168};
  Mat43 ACt = A * Transpose(C);
  for (int i = 0; i < ACt.Size(); ++i) {
    EXPECT_NEAR(ACt[i], ACt_expected[i], EPS);
  }

}
