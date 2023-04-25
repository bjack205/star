//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Mat3.hpp"

extern "C" {
#include "star/matrix3.h"
}

namespace star {

/*-------------------------------------
 * Static Methods
 *-----------------------------------*/

Mat3 Mat3::Zero() {
  Mat3 mat;
  star_SetZero33(mat.data());
  return mat;
}

Mat3::Mat3(sfloat x00, sfloat x10, sfloat x20, sfloat x01, sfloat x11, sfloat x21,
           sfloat x02, sfloat x12, sfloat x22)
    : data_{x00, x10, x20, x01, x11, x21, x02, x12, x22} {}

Mat3 Mat3::ByRows(sfloat x00, sfloat x01, sfloat x02, sfloat x10, sfloat x11, sfloat x12,
                  sfloat x20, sfloat x21, sfloat x22) {
  Mat3 mat;
  mat.data_[0] = x00;
  mat.data_[1] = x10;
  mat.data_[2] = x20;
  mat.data_[3] = x01;
  mat.data_[4] = x11;
  mat.data_[5] = x21;
  mat.data_[6] = x02;
  mat.data_[7] = x12;
  mat.data_[8] = x22;
  return mat;
}

/*-------------------------------------
 * Getters
 *-----------------------------------*/
Vec3 Mat3::GetRow(int row) const { return {data_[row], data_[row + 3], data_[row + 6]}; }

Vec3 Mat3::GetCol(int col) const {
  return {data_[col * 3], data_[col * 3 + 1], data_[col * 3 + 2]};
}

Vec3 Mat3::GetDiagonal() const { return {data_[0], data_[4], data_[8]}; }

/*-------------------------------------
 * Setters
 *-----------------------------------*/
void Mat3::SetRow(int row, const Vec3& v) {
  data_[row] = v[0];
  data_[row + 3] = v[1];
  data_[row + 6] = v[2];
}

void Mat3::SetCol(int col, const Vec3& v) {
  data_[col * 3] = v[0];
  data_[col * 3 + 1] = v[1];
  data_[col * 3 + 2] = v[2];
}
void Mat3::SetZero() { star_SetZero33(data_); }
void Mat3::SetIdentity() { star_SetIdentity33(data_, 1); }
void Mat3::SetConst(sfloat value) { star_SetConst33(data_, value); }

void Mat3::SetDiagonal(sfloat value) { star_SetIdentity33(data_, value); }

void Mat3::SetDiagonal(const Mat3& m) {
  Vec3 diag = m.GetDiagonal();
  star_SetDiagonal33(data_, diag.data());
}

void Mat3::SetDiagonal(sfloat x, sfloat y, sfloat z) {
  sfloat diag[3] = {x, y, z};
  star_SetDiagonal33(data_, diag);
}

Mat3 Mat3::Identity() {
  Mat3 mat;
  star_SetIdentity33(mat.data(), 1);
  return mat;
}

Mat3 Mat3::Const(sfloat value) {
  Mat3 mat;
  star_SetConst33(mat.data(), value);
  return mat;
}

Mat3 Mat3::Diagonal(sfloat value) {
  Mat3 mat;
  star_SetIdentity33(mat.data(), value);
  return mat;
}

Mat3 Mat3::Diagonal(const Mat3& m) {
  Mat3 mat;
  Vec3 diag = m.GetDiagonal();
  star_SetDiagonal33(mat.data(), diag.data());
  return mat;
}

Mat3 Mat3::Diagonal(sfloat x, sfloat y, sfloat z) {
  Mat3 mat;
  sfloat diag[3] = {x, y, z};
  star_SetDiagonal33(mat.data(), diag);
  return mat;
}

Mat3 Mat3::Transpose() const {
  Mat3 mat;
  star_Transpose33(mat.data(), data_);
  return mat;
}

Mat3& Mat3::TransposeInPlace() {
  star_TransposeInPlace33(data_);
  return *this;
}

}  // namespace star
