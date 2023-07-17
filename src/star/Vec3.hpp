//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <cstddef>

#include "star/typedefs.h"

namespace star {

class Vec3 {
 public:
  /*---------------------------------*/
  /* Constructors                    */
  /*---------------------------------*/
  Vec3() = default;
  Vec3(sfloat x, sfloat y, sfloat z) : x(x), y(y), z(z) {}

  template <class Vector>
  explicit Vec3(Vector v) : x(v[0]), y(v[1]), z(v[2]) {}

  /*---------------------------------*/
  /* Static Methods                  */
  /*---------------------------------*/
  static Vec3 Zero() { return {0, 0, 0}; }
  static Vec3 Const(sfloat value) { return {value, value, value}; }
  static Vec3 XAxis() { return {1, 0, 0}; }
  static Vec3 YAxis() { return {0, 1, 0}; }
  static Vec3 ZAxis() { return {0, 0, 1}; }

  /*---------------------------------*/
  /* Setters                         */
  /*---------------------------------*/
  void SetConst(sfloat value);
  void SetZero();

  /*---------------------------------*/
  /* Norms and Related               */
  /*---------------------------------*/
  sfloat Norm() const;
  sfloat NormSquared() const;
  sfloat InfNorm() const;
  sfloat OneNorm() const;
  Vec3 Normalize() const;
  Vec3& NormalizeInPlace();
  sfloat Dot(const Vec3& y) const;
  sfloat NormedDifference(const Vec3& other) const;

  Vec3 Cross(const Vec3& y) const;

  /*---------------------------------*/
  /* Element-wise operations         */
  /*---------------------------------*/
  Vec3 Add(const Vec3& y) const;
  Vec3 Sub(const Vec3& y) const;
  Vec3 Mul(const Vec3& y) const;
  Vec3 Div(const Vec3& y) const;

  Vec3& AddInPlace(const Vec3& y);
  Vec3& SubInPlace(const Vec3& y);
  Vec3& MulInPlace(const Vec3& y);
  Vec3& DivInPlace(const Vec3& y);

  Vec3 UnaryMap(sfloat (*function)(sfloat)) const;
  Vec3 BinaryMap(const Vec3& y, sfloat (*function)(sfloat, sfloat)) const;

  Vec3 operator+(const Vec3& rhs) const { return this->Add(rhs); }
  Vec3 operator-(const Vec3& rhs) const { return this->Sub(rhs); }
  Vec3 operator*(const Vec3& rhs) const { return this->Mul(rhs); }
  Vec3 operator/(const Vec3& rhs) const { return this->Div(rhs); }

  Vec3& operator+=(const Vec3& rhs) { return this->AddInPlace(rhs); };
  Vec3& operator-=(const Vec3& rhs) { return this->SubInPlace(rhs); };
  Vec3& operator*=(const Vec3& rhs) { return this->MulInPlace(rhs); };
  Vec3& operator/=(const Vec3& rhs) { return this->DivInPlace(rhs); };

  Vec3 operator+(sfloat alpha) const { return this->Add({alpha, alpha, alpha}); }
  Vec3 operator-(sfloat alpha) const { return this->Sub({alpha, alpha, alpha}); }
  Vec3 operator*(sfloat alpha) const { return this->Mul({alpha, alpha, alpha}); }
  Vec3 operator/(sfloat alpha) const { return this->Div({alpha, alpha, alpha}); }

  template <class T>
  Vec3& operator+=(T alpha) {
    x += alpha;
    y += alpha;
    z += alpha;
    return *this;
  }
  template <class T>
  Vec3& operator-=(T alpha) {
    x -= alpha;
    y -= alpha;
    z -= alpha;
    return *this;
  }
  template <class T>
  Vec3& operator*=(T alpha) {
    x *= alpha;
    y *= alpha;
    z *= alpha;
    return *this;
  }
  template <class T>
  Vec3& operator/=(T alpha) {
    x /= alpha;
    y /= alpha;
    z /= alpha;
    return *this;
  }

  /*---------------------------------*/
  /* Data Access                     */
  /*---------------------------------*/
  sfloat* data() { return &x; }
  const sfloat* data() const { return &x; }
  sfloat& operator[](size_t index) { return (&x)[index]; }
  const sfloat& operator[](size_t index) const { return (&x)[index]; }

  sfloat x;
  sfloat y;
  sfloat z;
};

static inline Vec3 operator+(sfloat alpha, const Vec3& rhs) { return rhs + alpha; }
static inline Vec3 operator-(sfloat alpha, const Vec3& rhs) { return rhs - alpha; }
static inline Vec3 operator*(sfloat alpha, const Vec3& rhs) { return rhs * alpha; }
static inline Vec3 operator/(sfloat alpha, const Vec3& rhs) { return rhs / alpha; }


}  // namespace star