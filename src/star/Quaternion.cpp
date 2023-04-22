//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Quaternion.hpp"

extern "C" {
#include "quaternion.h"
}

namespace star {

Quaternion::Quaternion(Vec4 v) : Vec4(std::move(v)) {}

/*---------------------------------*/
/* Static Methods                  */
/*---------------------------------*/
Quaternion Quaternion::Expm(sfloat x, sfloat y, sfloat z) {
  Quaternion q;
  sfloat phi[3] = {x, y, z};
  star_QuatExpm(q.data(), phi);
  return q;
}

Quaternion Quaternion::Expm(const Vec3& v) { return Expm(v[0], v[1], v[2]); }

Quaternion Quaternion::FromAxisAngle(sfloat angle, sfloat x, sfloat y, sfloat z) {
  Quaternion q;
  sfloat phi[3] = {angle * x, angle * y, angle * z};
  star_QuatExpm(q.data(), phi);
  return q;
}

Quaternion Quaternion::RotX(sfloat angle) {
  Quaternion q;
  star_QuatRotX(q.data(), angle);
  return q;
}

Quaternion Quaternion::RotY(sfloat angle) {
  Quaternion q;
  star_QuatRotY(q.data(), angle);
  return q;
}

Quaternion Quaternion::RotZ(sfloat angle) {
  Quaternion q;
  star_QuatRotZ(q.data(), angle);
  return q;
}

/*---------------------------------*/
/* Scalar Values                   */
/*---------------------------------*/
sfloat Quaternion::VecNorm() const { return star_QuatVecNorm(data()); }
sfloat Quaternion::VecNormSquared() const { return star_QuatVecNormSquared(data()); }
sfloat Quaternion::AngleBetween(const Quaternion rhs) const {
  return star_QuatAngleBetween(data(), rhs.data());
}

/*---------------------------------*/
/* Mathematical operators          */
/*---------------------------------*/
Quaternion Quaternion::Exp() const {
  Quaternion q;
  star_QuatExp(q.data(), data());
  return q;
}

Quaternion Quaternion::Log() const {
  Quaternion q;
  star_QuatLog(q.data(), data());
  return q;
}

Quaternion Quaternion::Flip() const {
  Quaternion q;
  star_QuatFlip(q.data(), data());
  return q;
}

Quaternion Quaternion::Conjugate() const {
  Quaternion q;
  star_QuatConjugate(q.data(), data());
  return q;
}

Quaternion Quaternion::Inverse() const {
  Quaternion q;
  star_QuatInverse(q.data(), data());
  return q;
}

Quaternion Quaternion::Compose(const Quaternion rhs) const {
  Quaternion q;
  star_QuatCompose(q.data(), data(), rhs.data());
  return q;
}

Quaternion Quaternion::ComposeLeft(const Quaternion lhs) const {
  Quaternion q;
  star_QuatComposeLeft(q.data(), lhs.data(), data());
  return q;
}

/*---------------------------------*/
/* Vector operations               */
/*---------------------------------*/

Vec3 Quaternion::RotateActive(const Vec3& v) const {
  Vec3 out;
  star_QuatRotateActive(out.data(), data(), v.data());
  return out;
}

Vec3 Quaternion::RotatePassive(const Vec3& v) const {
  Vec3 out;
  star_QuatRotatePassive(out.data(), data(), v.data());
  return out;
}

Quaternion Quaternion::ComposePure(const Vec3& v) const {
  Quaternion q;
  star_QuatComposePure(q.data(), data(), v.data());
  return q;
}

/*---------------------------------*/
/* Comparison                      */
/*---------------------------------*/

bool Quaternion::IsApprox(const Quaternion& rhs, sfloat tol) const {
  return star_QuatAngleBetween(data(), rhs.data()) < tol;
}

}  // namespace star
