//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#include "Vec4.hpp"

#include "star/vector4.h"

namespace star {

void Vec4::SetConst(sfloat value) { star_SetConst4(data(), value); }
void Vec4::SetZero() { star_SetZero4(data()); }
sfloat Vec4::Norm() const { return star_Norm4(data()); }
sfloat Vec4::NormSquared() const { return star_NormSquared4(data()); }
sfloat Vec4::InfNorm() const { return star_InfNorm4(data()); }
sfloat Vec4::OneNorm() const { return star_OneNorm4(data()); }

Vec4 Vec4::Normalize() const {
  Vec4 out(w, x, y, z);
  star_Normalize4(out.data(), data());
  return out;
}
Vec4& Vec4::NormalizeInPlace() {
  star_Normalize4(data(), data());
  return *this;
}

sfloat Vec4::Dot(const Vec4& y) const { return star_Dot4(data(), y.data()); }

sfloat Vec4::NormedDifference(const Vec4& other) const {
  return (*this - other).Norm();
}


Vec4 Vec4::Add(const Vec4& y) const {
  Vec4 out;
  star_Add4(out.data(), data(), y.data());
  return out;
}

Vec4 Vec4::Sub(const Vec4& y) const {
  Vec4 out;
  star_Sub4(out.data(), data(), y.data());
  return out;
}

Vec4 Vec4::Mul(const Vec4& y) const {
  Vec4 out;
  star_Mul4(out.data(), data(), y.data());
  return out;
}

Vec4 Vec4::Div(const Vec4& y) const {
  Vec4 out;
  star_Div4(out.data(), data(), y.data());
  return out;
}

Vec4& Vec4::AddInPlace(const Vec4& y) {
  star_Add4(data(), data(), y.data());
  return *this;
}

Vec4& Vec4::SubInPlace(const Vec4& y) {
  star_Sub4(data(), data(), y.data());
  return *this;
}

Vec4& Vec4::MulInPlace(const Vec4& y) {
  star_Mul4(data(), data(), y.data());
  return *this;
}

Vec4& Vec4::DivInPlace(const Vec4& y) {
  star_Div4(data(), data(), y.data());
  return *this;
}

}  // namespace star