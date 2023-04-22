//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "star/Transpose.hpp"
#include "star/Mat3.hpp"
#include "star/Vec3.hpp"

namespace star {

Mat3 Multiply(const Mat3& A, const Mat3& B);
Mat3 Multiply(const Transpose<Mat3>& At, const Mat3& B);
Mat3 Multiply(const Mat3& A, const Transpose<Mat3>& Bt);

Vec3 Multiply(const Mat3& A, const Vec3& B);
Vec3 Multiply(const Transpose<Mat3>& At, const Vec3& B);

void MultiplyInPlace(Mat3& C, const Mat3& A, const Mat3& B);

template <class Atype, class Btype>
auto operator*(const Atype& A, const Btype& B) {
  return Multiply(A, B);
}


}  // namespace star