//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "typedefs.h"
#include "Vec3.hpp"

namespace star {

class Mat3 {
 public:
  // Size information
  static constexpr int kRows = 3;
  static constexpr int kCols = 3;
  static constexpr int kSize = 9;
  constexpr int Rows() const { return kRows; }
  constexpr int Cols() const { return kCols; }
  constexpr int Size() const { return kSize; }

  // Constructors
  Mat3() = default;
  Mat3(sfloat x00, sfloat x10, sfloat x20, sfloat x01, sfloat x11, sfloat x21,
       sfloat x02, sfloat x12, sfloat x22);
  Mat3(const Mat3& other) = default;
  Mat3& operator=(const Mat3& other) = default;
  Mat3(Mat3&& other) = default;
  Mat3& operator=(Mat3&& other) = default;
  ~Mat3() = default;

  template <class Vector>
  explicit Mat3(Vector v) {
    for (int i = 0; i < kSize; ++i) {
      data_[i] = v[i];
    }
  }

  // Static methods
  static Mat3 ByRows(sfloat x00, sfloat x01, sfloat x02, sfloat x10, sfloat x11, sfloat x12,
                     sfloat x20, sfloat x21, sfloat x22);
  static Mat3 Zero();
  static Mat3 Identity();
  static Mat3 Identity(sfloat value);
  static Mat3 Const(sfloat value);
  static Mat3 Diagonal(sfloat value);
  static Mat3 Diagonal(const Mat3& m);
  static Mat3 Diagonal(sfloat x, sfloat y, sfloat z);

  template <class Vector>
  static Mat3 Diagonal(Vector v) {
    return Mat3::Diagonal(v[0], v[1], v[2]);
  }

  // Setters
  void SetZero();
  void SetIdentity();
  void SetIdentity(sfloat value);
  void SetConst(sfloat value);
  void SetDiagonal(sfloat value);
  void SetDiagonal(const Mat3& m);
  void SetDiagonal(sfloat x, sfloat y, sfloat z);

  // Transpose
  Mat3 Transpose() const;
  Mat3& TransposeInPlace();

//  // Multiplication
//  template <class MatOrVec>
//  auto Mul(const MatOrVec& B) const {
//    return Multiply(*this, B);
//  }
//
//  template <class MatOrVec>
//  auto operator*(const MatOrVec& B) const {
//    return Multiply(*this, B);
//  }

  // Data Access
  sfloat& operator[](int i) { return data_[i]; }
  const sfloat& operator[](int i) const { return data_[i]; }
  sfloat& operator()(int i, int j) { return data_[i * kCols + j]; }
  const sfloat& operator()(int i, int j) const { return data_[i * kCols + j]; }
  sfloat* data() { return data_; }
  const sfloat* data() const { return data_; }


 private:
  sfloat data_[kSize];
};

}  // namespace star
