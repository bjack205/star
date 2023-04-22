//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Mat3.hpp"

extern "C" {
#include "star/matrix3.h"
}

namespace star {

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

}  // namespace star
