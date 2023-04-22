//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>

#include "star/Vec3.hpp"
#include "star/Vec4.hpp"
#include "star/typedefs.h"

namespace star {

class Matrix3;
class Matrix4;
class Matrix43;
class Matrix34;

class Quaternion : public Vec4 {
 public:
  using Vec4::Vec4;

  /*---------------------------------*/
  /* Constructors                    */
  /*---------------------------------*/
  Quaternion() : Vec4(1,0,0,0) {}
  Quaternion(sfloat w, sfloat i, sfloat j, sfloat k) : Vec4(w, i, j, k) {}
  Quaternion(sfloat w, const Vec3& v) : Vec4(w, v[0], v[1], v[2]) {}
  Quaternion(Vec4 v);  // NOLINT: Allow implicit conversion

  template <class Vector>
  explicit Quaternion(const Vector& v) : Vec4(v) {}

  /*---------------------------------*/
  /* Static Methods                  */
  /*---------------------------------*/
  static Quaternion Identity() { return {1, 0, 0, 0}; };
  static Quaternion Pure(sfloat x, sfloat y, sfloat z) { return {0, x, y, z}; };
  static Quaternion Pure(const Vec3& v) { return {0, v[0], v[1], v[2]}; };
  static Quaternion Expm(sfloat x, sfloat y, sfloat z);
  static Quaternion Expm(const Vec3& v);
  static Quaternion FromAxisAngle(sfloat angle, sfloat x, sfloat y, sfloat z);
  static Quaternion RotX(sfloat angle);
  static Quaternion RotY(sfloat angle);
  static Quaternion RotZ(sfloat angle);

  /*---------------------------------*/
  /* Scalar Values                   */
  /*---------------------------------*/
  sfloat VecNorm() const;
  sfloat VecNormSquared() const;
  sfloat AngleBetween(const Quaternion rhs) const;

  /*---------------------------------*/
  /* Mathematical operators          */
  /*---------------------------------*/
  Quaternion Exp() const;
  Quaternion Log() const;
  Quaternion Flip() const;
  Quaternion Conjugate() const;
  Quaternion Inverse() const;
  Quaternion Compose(const Quaternion rhs) const;
  Quaternion ComposeLeft(const Quaternion lhs) const;

  /*---------------------------------*/
  /* Vector operations               */
  /*---------------------------------*/
  Vec3 Vec() const { return { i, j, k }; };
  Vec3 RotateActive(const Vec3& v) const;
  Vec3 RotatePassive(const Vec3& v) const;
  Quaternion ComposePure(const Vec3& v) const;

  /*---------------------------------*/
  /* Jacobians                       */
  /*---------------------------------*/
  Matrix34 RotActiveJacobian(const Vec3&) const;
  Matrix34 RotPassiveJacobian(const Vec3&) const;
  Matrix43 AttitudeJacobian() const;

  /*---------------------------------*/
  /* Matrices                       */
  /*---------------------------------*/
  Matrix4 L() const;
  Matrix4 R() const;

  /*---------------------------------*/
  /* Comparison                      */
  /*---------------------------------*/
  bool IsApprox(const Quaternion& rhs, sfloat tol = 1e-6) const;
};

}  // namespace star
