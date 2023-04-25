//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "matrix_multiplication.hpp"

extern "C" {
#include "star/matrix3.h"
#include "star/matrix4.h"
#include "star/matrix43.h"
}

namespace star {

/*-------------------------------------
 * 3x3 Matrices
 *-----------------------------------*/
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

Vec3 Multiply(const Transpose<Mat3>& At, const Vec3& x) {
  Vec3 y;
  star_TransposedVecMul33(y.data(), At.data(), x.data());
  return y;
}

void MultiplyInPlace(Mat3& C, const Mat3& A, const Mat3& B) {
  star_MatMul33(C.data(), A.data(), B.data());
}

void MultiplyInPlace(Mat3& C, const Transpose<Mat3>& At, const Mat3& B) {
  star_TransposedMatMul33(C.data(), At.data(), B.data());
}

void MultiplyInPlace(Mat3& C, const Mat3& A, const Transpose<Mat3>& B) {
  star_MatMulTransposed33(C.data(), A.data(), B.data());
}

/*-------------------------------------
 * 4x4 Matrices
 *-----------------------------------*/

Mat4 Multiply(const Mat4& A, const Mat4& B) {
  Mat4 C;
  star_MatMul44(C.data(), A.data(), B.data());
  return C;
}

Mat4 Multiply(const Transpose<Mat4>& At, const Mat4& B) {
  Mat4 C;
  star_TransposedMatMul44(C.data(), At.data(), B.data());
  return C;
}

Mat4 Multiply(const Mat4& A, const Transpose<Mat4>& Bt) {
  Mat4 C;
  star_MatMulTransposed44(C.data(), A.data(), Bt.data());
  return C;
}

Vec4 Multiply(const Mat4& A, const Vec4& x) {
  Vec4 y;
  star_VecMul44(y.data(), A.data(), x.data());
  return y;
}

Vec4 Multiply(const Transpose<Mat4>& At, const Vec4& x) {
  Vec4 y;
  star_TransposedVecMul44(y.data(), At.data(), x.data());
  return y;
}

void MultiplyInPlace(Mat4& C, const Mat4& A, const Mat4& B) {
  star_MatMul44(C.data(), A.data(), B.data());
}

void MultiplyInPlace(Mat4& C, const Transpose<Mat4>& At, const Mat4& B) {
  star_TransposedMatMul44(C.data(), At.data(), B.data());
}

void MultiplyInPlace(Mat4& C, const Mat4& A, const Transpose<Mat4>& B) {
  star_MatMulTransposed44(C.data(), A.data(), B.data());
}

/*-------------------------------------
 * 4x3 Matrices
 *-----------------------------------*/

Mat43 Multiply(const Mat43& A, const Mat3& B) {
  Mat43 C;
  star_MatMul433(C.data(), A.data(), B.data());
  return C;
}

Mat43 Multiply(const Mat4& A, const Mat43& B) {
  Mat43 C;
  star_MatMul443(C.data(), A.data(), B.data());
  return C;
}

Mat43 Multiply(const Mat43& A, const Transpose<Mat3>& B) {
  Mat43 C;
  star_MatMulTransposed433(C.data(), A.data(), B.data());
  return C;
}

Mat43 Multiply(const Transpose<Mat4>& A, const Mat43& B) {
  Mat43 C;
  star_TransposedMatMul443(C.data(), A.data(), B.data());
  return C;
}

Vec4 Multiply(const Mat43& A, const Vec3& x) {
  Vec4 y;
  star_VecMul43(y.data(), A.data(), x.data());
  return y;
}

Vec3 Multiply(const Transpose<Mat43>& A, const Vec4& x) {
  Vec3 y;
  star_TransposedVecMul43(y.data(), A.data(), x.data());
  return y;
}

void MultiplyInPlace(Mat43& C, const Mat43& A, const Mat3& B) {
  star_MatMul433(C.data(), A.data(), B.data());
}

void MultiplyInPlace(Mat43& C, const Mat4& A, const Mat43& B) {
  star_MatMul443(C.data(), A.data(), B.data());
}

void MultiplyInPlace(Mat43& C, const Transpose<Mat4>& A, const Mat43& B) {
  star_TransposedMatMul443(C.data(), A.data(), B.data());
}

void MultiplyInPlace(Vec4& y, const Mat43& A, const Vec3& x) {
  star_VecMul43(y.data(), A.data(), x.data());
}

void MultiplyInPlace(Vec3& y, const Transpose<Mat43>& A, const Vec4& x) {
  star_TransposedVecMul43(y.data(), A.data(), x.data());
}

//void MultiplyInPlace(Transpose<Mat43>& C, const Mat3& A, const Transpose<Mat43>& B) {
//  star_MatMul334(C.data(), A.data(), B.data());
//}
//
//void MultiplyInPlace(Transpose<Mat43>& C, const Transpose<Mat43>& A, const Mat4& B) {
//  star_MatMul344(C.data(), A.data(), B.data());
//}

}  // namespace star
