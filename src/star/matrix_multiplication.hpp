//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "star/Mat3.hpp"
#include "star/Mat4.hpp"
#include "star/Mat43.hpp"
#include "star/Transpose.hpp"
#include "star/Vec3.hpp"
#include "star/Vec4.hpp"

namespace star {

/*
 * Generic multiplication operator overload
 */
template <class Atype, class Btype>
auto operator*(const Atype& A, const Btype& B) {
  return Multiply(A, B);
}

/*-------------------------------------
 * 3x3 Matrices
 *-----------------------------------*/
Mat3 Multiply(const Mat3& A, const Mat3& B);
Mat3 Multiply(const Transpose<Mat3>& At, const Mat3& B);
Mat3 Multiply(const Mat3& A, const Transpose<Mat3>& Bt);

Vec3 Multiply(const Mat3& A, const Vec3& x);
Vec3 Multiply(const Transpose<Mat3>& At, const Vec3& x);

void MultiplyInPlace(Mat3& C, const Mat3& A, const Mat3& B);
void MultiplyInPlace(Mat3& C, const Transpose<Mat3>& At, const Mat3& B);
void MultiplyInPlace(Mat3& C, const Mat3& A, const Transpose<Mat3>& B);

/*-------------------------------------
 * 4x4 Matrices
 *-----------------------------------*/
Mat4 Multiply(const Mat4& A, const Mat4& B);
Mat4 Multiply(const Transpose<Mat4>& At, const Mat4& B);
Mat4 Multiply(const Mat4& A, const Transpose<Mat4>& Bt);

Vec4 Multiply(const Mat4& A, const Vec4& x);
Vec4 Multiply(const Transpose<Mat4>& At, const Vec4& x);

void MultiplyInPlace(Mat4& C, const Mat4& A, const Mat4& B);
void MultiplyInPlace(Mat4& C, const Transpose<Mat4>& At, const Mat4& B);
void MultiplyInPlace(Mat4& C, const Mat4& A, const Transpose<Mat4>& B);

/*-------------------------------------
 * 4x3 Matrices
 *-----------------------------------*/
Mat43 Multiply(const Mat43& A, const Mat3& B);
Mat43 Multiply(const Mat4& A, const Mat43& B);

Mat43 Multiply(const Mat43& A, const Transpose<Mat3>& B);
Mat43 Multiply(const Transpose<Mat4>& A, const Mat43& B);

Vec4 Multiply(const Mat43& A, const Vec3& x);
Vec3 Multiply(const Transpose<Mat43>& A, const Vec4& x);

void MultiplyInPlace(Mat43& C, const Mat43& A, const Mat3& B);
void MultiplyInPlace(Mat43& C, const Mat4& A, const Mat43& B);

void MultiplyInPlace(Mat43& C, const Mat43& A, const Transpose<Mat3>& B);
void MultiplyInPlace(Mat43& C, const Transpose<Mat4>& A, const Mat43& B);

void MultiplyInPlace(Transpose<Mat43>& C, const Mat3& A, const Transpose<Mat43>& B);
void MultiplyInPlace(Transpose<Mat43>& C, const Transpose<Mat43>& A, const Mat4& B);

void MultiplyInPlace(Vec4& y, const Mat43& A, const Vec3& x);
void MultiplyInPlace(Vec3& y, const Transpose<Mat43>& A, const Vec4& x);

}  // namespace star