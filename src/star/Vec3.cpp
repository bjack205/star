//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Vec3.hpp"

#include "star/vector3.h"

namespace star {

sfloat Vec3::Norm() const { return star_Norm3(&x); }
sfloat Vec3::NormSquared() const { return star_NormSquared3(&x); }
sfloat Vec3::InfNorm() const { return star_InfNorm3(&x); }
sfloat Vec3::OneNorm() const { return star_OneNorm3(&x); }

Vec3 Vec3::Normalize() const {
  Vec3 out(x, y, z);
  star_Normalize3(out.data(), data());
  return out;
}

Vec3& Vec3::NormalizeInPlace() {
  star_Normalize3(data(), data());
  return *this;
}

sfloat Vec3::Dot(const Vec3& y) const { return star_Dot3(data(), y.data()); }

Vec3 Vec3::Add(const Vec3& y) const {
  Vec3 out;
  star_Add3(out.data(), data(), y.data());
  return out;
}

Vec3 Vec3::Sub(const Vec3& y) const {
  Vec3 out;
  star_Sub3(out.data(), data(), y.data());
  return out;
}

Vec3 Vec3::Mul(const Vec3& y) const {
  Vec3 out;
  star_Mul3(out.data(), data(), y.data());
  return out;
}

Vec3 Vec3::Div(const Vec3& y) const {
  Vec3 out;
  star_Div3(out.data(), data(), y.data());
  return out;
}

Vec3& Vec3::AddInPlace(const Vec3& y) {
  star_Add3(data(), data(), y.data());
  return *this;
}

Vec3& Vec3::SubInPlace(const Vec3& y) {
  star_Sub3(data(), data(), y.data());
  return *this;
}

Vec3& Vec3::MulInPlace(const Vec3& y) {
  star_Mul3(data(), data(), y.data());
  return *this;
}

Vec3& Vec3::DivInPlace(const Vec3& y) {
  star_Div3(data(), data(), y.data());
  return *this;
}

Vec3 Vec3::UnaryMap(sfloat (*function)(sfloat)) const {
  Vec3 out;
  star_UnaryMap(out.data(), data(), function);
  return out;
}

Vec3 Vec3::BinaryMap(const star::Vec3& y, sfloat (*function)(sfloat, sfloat)) const {
  Vec3 out;
  star_BinaryMap(out.data(), data(), y.data(), function);
  return out;
}

void Vec3::SetConst(sfloat value) {
  star_SetConst(data(), value);
}

void Vec3::SetZero() {
  star_SetZero(data());
}

}  // namespace star