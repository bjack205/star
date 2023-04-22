//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "matrix_multiplication.hpp"

extern "C" {
#include "star/matrix3.h"
}

namespace star {

Mat3 Multiply(const Mat3& A, const Mat3& B) {
  Mat3 C;
  star_MatMul33(C.data(), A.data(), B.data());
  return C;
}

Mat3 Multiply(const Transpose<Mat3>& At, const Mat3& B) {
  Mat3 C;
  star_TransposedMatMul33(C.data(), At.data(), B.data());
  return C;
}

Mat3 Multiply(const Mat3& A, const Transpose<Mat3>& Bt) {
  Mat3 C;
  star_MatMulTransposed33(C.data(), A.data(), Bt.data());
  return C;
}

Vec3 Multiply(const Mat3& A, const Vec3& x) {
  Vec3 y;
  star_VecMul33(y.data(), A.data(), x.data());
  return y;
}

Vec3 Multiply(const Transpose<Mat3>& At, const Vec3& B) {
  Vec3 y;
  star_TransposedVecMul33(y.data(), At.data(), B.data());
  return y;
}

void MultiplyInPlace(Mat3& C, const Mat3& A, const Mat3& B) {
  star_MatMul33(C.data(), A.data(), B.data());
}

}  // namespace star
