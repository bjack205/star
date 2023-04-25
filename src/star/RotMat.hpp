//
// Created by Brian Jackson on 4/25/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <cmath>
#include "star/Mat3.hpp"
#include "star/typedefs.h"

namespace star {

constexpr bool Active = true;
constexpr bool Passive = false;

template <bool Sense>
class RotMat : public Mat3 {
 public:
  /*-------------------------------------
   * Constructors
   *-----------------------------------*/
  // NOTE: RotMat does NOT inherit the constructors of Mat3
  RotMat() : Mat3(Identity()) {}
  explicit RotMat(const Mat3& mat) : Mat3(mat) {}
  RotMat(sfloat R00, sfloat R01, sfloat R02, sfloat R10, sfloat R11, sfloat R12, sfloat R20,
         sfloat R21, sfloat R22)
      : Mat3(R00, R01, R02, R10, R11, R12, R20, R21, R22) {}

  /*-------------------------------------
   * Static Methods
   *-----------------------------------*/
  static RotMat ByRows(sfloat R00, sfloat R01, sfloat R02, sfloat R10, sfloat R11,
                       sfloat R12, sfloat R20, sfloat R21, sfloat R22) {
    return RotMat(R00, R01, R02, R10, R11, R12, R20, R21, R22);
  }

  static RotMat RotX(sfloat angle);
  static RotMat RotY(sfloat angle);
  static RotMat RotZ(sfloat angle);
  static RotMat FromAxisAngle(sfloat angle, sfloat x, sfloat y, sfloat z);
  static RotMat FromAxisAngle(sfloat angle, const Vec3& axis);
  static RotMat FromQuaternion(sfloat w, sfloat i, sfloat j, sfloat k);


  /*-------------------------------------
   * Linear Algebra
   *-----------------------------------*/
  RotMat Transpose() const { return RotMat(Mat3::Transpose()); }
  RotMat& TransposeInPlace() {
    Mat3::TransposeInPlace();
    return *this;
  }
};

/*-------------------------------------
 * Active Rotations
 *-----------------------------------*/
template <>
RotMat<Active> RotMat<Active>::RotX(sfloat angle) {
  sfloat c = std::cos(angle);
  sfloat s = std::sin(angle);
  return RotMat<Active>(1, 0, 0, 0, c, -s, 0, s, c);
}

template <>
RotMat<Active> RotMat<Active>::RotY(sfloat angle) {
  sfloat c = std::cos(angle);
  sfloat s = std::sin(angle);
  return RotMat<Active>(c, 0, s, 0, 1, 0, -s, 0, c);
}

template <>
RotMat<Active> RotMat<Active>::RotZ(sfloat angle) {
  sfloat c = std::cos(angle);
  sfloat s = std::sin(angle);
  return RotMat<Active>(c, -s, 0, s, c, 0, 0, 0, 1);
}

/*-------------------------------------
 * Passive Rotations
 *-----------------------------------*/
template <>
RotMat<Passive> RotMat<Passive>::RotX(sfloat angle) {
  sfloat c = std::cos(angle);
  sfloat s = std::sin(angle);
  return RotMat<Passive>(1, 0, 0, 0, c, s, 0, -s, c);
}

template <>
RotMat<Passive> RotMat<Passive>::RotY(sfloat angle) {
  sfloat c = std::cos(angle);
  sfloat s = std::sin(angle);
  return RotMat<Passive>(c, 0, -s, 0, 1, 0, s, 0, c);
}

template <>
RotMat<Passive> RotMat<Passive>::RotZ(sfloat angle) {
  sfloat c = std::cos(angle);
  sfloat s = std::sin(angle);
  return RotMat<Passive>(c, s, 0, -s, c, 0, 0, 0, 1);
}

}  // namespace star