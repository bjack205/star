//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "star/Vec3.hpp"
#include "star/Vec4.hpp"
#include "star/typedefs.h"

namespace star {

class Mat43 {
 public:
  // Size information
  static constexpr int kRows = 4;
  static constexpr int kCols = 3;
  static constexpr int kSize = 12;
  constexpr int Rows() const { return kRows; }
  constexpr int Cols() const { return kCols; }
  constexpr int Size() const { return kSize; }

  /*-------------------------------------
   * Constructors
   *-----------------------------------*/
  Mat43() = default;
  Mat43(sfloat x00, sfloat x10, sfloat x20, sfloat x30, sfloat x01, sfloat x11, sfloat x21,
        sfloat x31, sfloat x02, sfloat x12, sfloat x22, sfloat x32)
      : data_{x00, x10, x20, x30, x01, x11, x21, x31, x02, x12, x22, x32} {}

  template <class Vector>
  explicit Mat43(Vector v) {
    for (int i = 0; i < kSize; ++i) {
      data_[i] = v[i];
    }
  }

  /*-------------------------------------
   * Static Methods
   *-----------------------------------*/
  static Mat43 ByRows(sfloat x00, sfloat x01, sfloat x02, sfloat x10, sfloat x11,
                      sfloat x12, sfloat x20, sfloat x21, sfloat x22, sfloat x30,
                      sfloat x31, sfloat x32);
  static Mat43 Zero();
  static Mat43 Const(sfloat value);

  /*-------------------------------------
   * Getters
   *-----------------------------------*/
  Vec3 GetRow(int i) const { return {data_[i * 3], data_[i * 3 + 1], data_[i * 3 + 2]}; }
  Vec4 GetCol(int j) const { return {data_[j], data_[j + 3], data_[j + 6], data_[j + 9]}; }

  /*-------------------------------------
   * Setters
   *-----------------------------------*/
  void SetRow(int i, const Vec3& v);
  void SetCol(int j, const Vec4& v);
  void SetZero();
  void SetConst(sfloat value);

  /*-------------------------------------
   * Data Access
   * -----------------------------------*/
  sfloat& operator[](int k) { return data_[k]; }
  const sfloat& operator[](int k) const { return data_[k]; }
  sfloat& operator[](IndexPair ij) { return data_[std::get<0>(ij) + kRows * std::get<1>(ij)]; }
  const sfloat& operator[](IndexPair ij) const {
    return data_[std::get<0>(ij) + kRows * std::get<1>(ij)];
  }

  sfloat& operator()(int i, int j) { return data_[i + 4 * j]; }
  const sfloat& operator()(int i, int j) const { return data_[i + 4 * j]; }
  sfloat* data() { return data_; }
  const sfloat* data() const { return data_; }

 private:
  sfloat data_[kSize] = {0};
};

}  // namespace star