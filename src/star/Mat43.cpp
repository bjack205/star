//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Mat43.hpp"

extern "C" {
#include "star/matrix43.h"
}

namespace star {

Mat43 Mat43::ByRows(sfloat x00, sfloat x01, sfloat x02, sfloat x10, sfloat x11, sfloat x12,
                    sfloat x20, sfloat x21, sfloat x22, sfloat x30, sfloat x31,
                    sfloat x32) {
  Mat43 mat;
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
  return mat;
}

Mat43 Mat43::Zero() {
  Mat43 mat;
  star_SetZero43(mat.data());
  return mat;
}

Mat43 Mat43::Const(sfloat value) {
  Mat43 mat;
  star_SetConst43(mat.data(), value);
  return mat;
}

void Mat43::SetRow(int i, const Vec3& v) {
  data_[i] = v[0];
  data_[i + 4] = v[1];
  data_[i + 8] = v[2];
}

void Mat43::SetCol(int j, const Vec4& v) {
  data_[j * 4] = v[0];
  data_[j * 4 + 1] = v[1];
  data_[j * 4 + 2] = v[2];
  data_[j * 4 + 3] = v[3];
}

void Mat43::SetZero() { star_SetZero43(data()); }
void Mat43::SetConst(sfloat value) { star_SetConst43(data(), value); }

/*-------------------------------------
 * Setters
 *-----------------------------------*/

}  // namespace star