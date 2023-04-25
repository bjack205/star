//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Mat4.hpp"

extern "C" {
#include "star/matrix4.h"
}

namespace star {

/*-------------------------------------
 * Constructors
 * -----------------------------------*/
Mat4::Mat4(sfloat x00, sfloat x10, sfloat x20, sfloat x30, sfloat x01, sfloat x11,
           sfloat x21, sfloat x31, sfloat x02, sfloat x12, sfloat x22, sfloat x32,
           sfloat x03, sfloat x13, sfloat x23, sfloat x33)
    : data_{
          x00, x10, x20, x30, x01, x11, x21, x31, x02, x12, x22, x32, x03, x13, x23, x33,
      } {}

/*-------------------------------------
 * Static Methods
 * -----------------------------------*/

Mat4 Mat4::ByRows(sfloat x00, sfloat x01, sfloat x02, sfloat x03, sfloat x10, sfloat x11,
                  sfloat x12, sfloat x13, sfloat x20, sfloat x21, sfloat x22, sfloat x23,
                  sfloat x30, sfloat x31, sfloat x32, sfloat x33) {
  Mat4 mat;
  mat.data_[0] = x00;
  mat.data_[1] = x10;
  mat.data_[2] = x20;
  mat.data_[3] = x30;
  mat.data_[4] = x01;
  mat.data_[5] = x11;
  mat.data_[6] = x21;
  mat.data_[7] = x31;
  mat.data_[8] = x02;
  mat.data_[9] = x12;
  mat.data_[10] = x22;
  mat.data_[11] = x32;
  mat.data_[12] = x03;
  mat.data_[13] = x13;
  mat.data_[14] = x23;
  mat.data_[15] = x33;
  return mat;
}

Mat4 Mat4::Zero() {
  Mat4 mat;
  star_SetZero44(mat.data());
  return mat;
}

Mat4 Mat4::Identity() {
  Mat4 mat;
  star_SetIdentity44(mat.data(), 1);
  return mat;
}

Mat4 Mat4::Const(sfloat value) {
  Mat4 mat;
  star_SetConst44(mat.data(), value);
  return mat;
}

Mat4 Mat4::Diagonal(sfloat value) {
  Mat4 mat;
  star_SetIdentity44(mat.data(), value);
  return mat;
}

Mat4 Mat4::Diagonal(sfloat x, sfloat y, sfloat z, sfloat w) {
  Mat4 mat;
  sfloat diag[4] = {x, y, z, w};
  star_SetDiagonal44(mat.data(), diag);
  return mat;
}

Mat4 Mat4::Diagonal(const Mat4& m) { return Diagonal(m[0], m[5], m[10], m[15]); }

/*-------------------------------------
 * Getters
 *-----------------------------------*/
Vec4 Mat4::GetRow(int row) const {
  return {data_[row], data_[row + 4], data_[row + 8], data_[row + 12]};
}

Vec4 Mat4::GetCol(int col) const {
  return {data_[col * 4], data_[col * 4 + 1], data_[col * 4 + 2], data_[col * 4 + 3]};
}

Vec4 Mat4::GetDiagonal() const { return {data_[0], data_[5], data_[10], data_[15]}; }

/*-------------------------------------
 * Setters
 *-----------------------------------*/

void Mat4::SetRow(int row, const Vec4& v) {
  data_[row] = v[0];
  data_[row + 4] = v[1];
  data_[row + 8] = v[2];
  data_[row + 12] = v[3];
}

void Mat4::SetCol(int col, const Vec4& v) {
  data_[col * 4] = v[0];
  data_[col * 4 + 1] = v[1];
  data_[col * 4 + 2] = v[2];
  data_[col * 4 + 3] = v[3];
}

void Mat4::SetZero() { star_SetZero44(data_); }

void Mat4::SetIdentity() { star_SetIdentity44(data_, 1); }

void Mat4::SetConst(sfloat value) { star_SetConst44(data_, value); }

void Mat4::SetDiagonal(sfloat value) { star_SetIdentity44(data_, value); }

void Mat4::SetDiagonal(sfloat x, sfloat y, sfloat z, sfloat w) {
  sfloat diag[4] = {x, y, z, w};
  star_SetDiagonal44(data_, diag);
}

void Mat4::SetDiagonal(const Mat4& m) { SetDiagonal(m[0], m[5], m[10], m[15]); }

/*-------------------------------------
 * Linear Algebra
 *-----------------------------------*/

Mat4 Mat4::Transpose() const {
  Mat4 mat;
  star_Transpose44(mat.data(), data_);
  return mat;
}

Mat4& Mat4::TransposeInPlace() {
  star_TransposeInPlace44(data_);
  return *this;
}

}  // namespace star
