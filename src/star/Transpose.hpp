//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "Mat3.hpp"
#include "Vec3.hpp"
#include "typedefs.h"

namespace star {

template <class Mat>
class Transpose {
  // Static assert that Mat is a matrix
  static_assert(Mat::kCols > 1, "Must be a matrix");

 public:
  // Size information
  static constexpr int kRows = Mat::kCols;
  static constexpr int kCols = Mat::kRows;
  static constexpr int kSize = Mat::kSize;
  constexpr int Rows() const { return kRows; }
  constexpr int Cols() const { return kCols; }
  constexpr int Size() const { return kSize; }

  /*-------------------------------------
   * Constructors
   *-----------------------------------*/
  Transpose(Mat& mat) : mat_(mat) {};

  /*-------------------------------------
   * Getters
   *-----------------------------------*/
  auto GetRow(int row) const {
    return mat_.GetCol(row);
  }

  auto GetCol(int col) const {
    return mat_.GetRow(col);
  }

  /*-------------------------------------
   * Data Access
   *-----------------------------------*/
  sfloat& operator[](int k) { return mat_(k % 3, k / 3); }
  const sfloat& operator[](int k) const { return mat_(k % 3, k / 3); }
  sfloat& operator()(int i, int j) { return mat_(j, i); }
  const sfloat& operator()(int i, int j) const { return mat_(j, i); }
  sfloat* data() { return mat_.data(); }
  const sfloat* data() const { return mat_.data(); }

 private:
  Mat& mat_;
};


}  // namespace star