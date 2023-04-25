//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "star/Vec4.hpp"
#include "typedefs.h"

namespace star {

class Mat4 {
 public:
  // Size Info
  static constexpr int kRows = 4;
  static constexpr int kCols = 4;
  static constexpr int kSize = 16;
  constexpr int Rows() const { return kRows; }
  constexpr int Cols() const { return kCols; }
  constexpr int Size() const { return kSize; }

  /*-------------------------------------
   * Constructors
   *-----------------------------------*/
  Mat4() = default;
  Mat4(sfloat x00, sfloat x10, sfloat x20, sfloat x30, sfloat x01, sfloat x11, sfloat x21,
       sfloat x31, sfloat x02, sfloat x12, sfloat x22, sfloat x32, sfloat x03, sfloat x13,
       sfloat x23, sfloat x33);

  template <class Vector>
  explicit Mat4(Vector v) {
    for (int i = 0; i < kSize; ++i) {
      data_[i] = v[i];
    }
  }

  /*-------------------------------------
   * Static Methods
   *-----------------------------------*/
  static Mat4 ByRows(sfloat x00, sfloat x01, sfloat x02, sfloat x03, sfloat x10, sfloat x11,
                     sfloat x12, sfloat x13, sfloat x20, sfloat x21, sfloat x22, sfloat x23,
                     sfloat x30, sfloat x31, sfloat x32, sfloat x33);
  static Mat4 Zero();
  static Mat4 Identity();
  static Mat4 Const(sfloat value);
  static Mat4 Diagonal(sfloat value);
  static Mat4 Diagonal(sfloat x, sfloat y, sfloat z, sfloat w);
  static Mat4 Diagonal(const Mat4& m);

  template <class Vector>
  static Mat4 Diagonal(Vector v) {
    return Mat4::Diagonal(v[0], v[1], v[2], v[3]);
  }

  /*-------------------------------------
   * Getters
   *-----------------------------------*/
  Vec4 GetRow(int row) const;
  Vec4 GetCol(int col) const;
  Vec4 GetDiagonal() const;

  /*-------------------------------------
   * Setters
   *-----------------------------------*/
  void SetRow(int row, const Vec4& v);
  void SetCol(int col, const Vec4& v);
  void SetZero();
  void SetIdentity();
  void SetConst(sfloat value);
  void SetDiagonal(sfloat value);
  void SetDiagonal(const Mat4& m);
  void SetDiagonal(sfloat x, sfloat y, sfloat z, sfloat w);

  template <class Vector>
  void SetDiagonal(Vector v) {
    SetDiagonal(v[0], v[1], v[2], v[3]);
  }

  /*-------------------------------------
   * Linear Algebra
   *-----------------------------------*/
  // TODO: Add linear algebra
  Mat4 Transpose() const;
  Mat4& TransposeInPlace();

  /*-------------------------------------
   * Data Access
   *-----------------------------------*/
  sfloat& operator[](int index) { return data_[index]; }
  const sfloat& operator[](int index) const { return data_[index]; }
  sfloat& operator[](IndexPair ij) { return data_[std::get<0>(ij) + kRows * std::get<1>(ij)]; }
  const sfloat& operator[](IndexPair ij) const {
    return data_[std::get<0>(ij) + kRows * std::get<1>(ij)];
  }

  sfloat& operator()(int row, int col) { return data_[row + col * kRows]; }
  const sfloat& operator()(int row, int col) const { return data_[row + col * kRows]; }
  sfloat* data() { return data_; }
  const sfloat* data() const { return data_; }

 private:
  sfloat data_[kSize];
};

}  // namespace star