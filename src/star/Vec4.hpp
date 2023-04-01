//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <cstddef>

#include "star/typedefs.h"

namespace star {

class Vec4 {
 public:
  /*---------------------------------*/
  /* Constructors                    */
  /*---------------------------------*/
  Vec4() = default;
  Vec4(sfloat w, sfloat x, sfloat y, sfloat z) : w(w),  x(x), y(y), z(z) {}

  template <class Vector>
  Vec4(Vector v) : w(v[0]), x(v[1]), y(v[2]), z(v[3]) {}

  /*---------------------------------*/
  /* Static Methods                  */
  /*---------------------------------*/
  static Vec4 Zero() { return {0, 0, 0, 0}; }
  static Vec4 Const(sfloat value) { return {value, value, value, value}; }

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
  Vec4 Normalize() const;
  Vec4& NormalizeInPlace();
  sfloat Dot(const Vec4& y) const;

  /*---------------------------------*/
  /* Element-wise operations         */
  /*---------------------------------*/
  Vec4 Add(const Vec4& y) const;
  Vec4 Sub(const Vec4& y) const;
  Vec4 Mul(const Vec4& y) const;
  Vec4 Div(const Vec4& y) const;

  Vec4& AddInPlace(const Vec4& y);
  Vec4& SubInPlace(const Vec4& y);
  Vec4& MulInPlace(const Vec4& y);
  Vec4& DivInPlace(const Vec4& y);

  Vec4 UnaryMap(sfloat (*function)(sfloat)) const;
  Vec4 BinaryMap(const Vec4& y, sfloat (*function)(sfloat, sfloat)) const;


  Vec4 operator+(const Vec4& rhs) const { return this->Add(rhs); }
  Vec4 operator-(const Vec4& rhs) const { return this->Sub(rhs); }
  Vec4 operator*(const Vec4& rhs) const { return this->Mul(rhs); }
  Vec4 operator/(const Vec4& rhs) const { return this->Div(rhs); }

  Vec4& operator+=(const Vec4& rhs) { return this->AddInPlace(rhs); };
  Vec4& operator-=(const Vec4& rhs) { return this->SubInPlace(rhs); };
  Vec4& operator*=(const Vec4& rhs) { return this->MulInPlace(rhs); };
  Vec4& operator/=(const Vec4& rhs) { return this->DivInPlace(rhs); };

  /*---------------------------------*/
  /* Data Access                     */
  /*---------------------------------*/
  sfloat *data() { return &w; }
  const sfloat *data() const { return &w; }
  sfloat& operator[](size_t index) { return (&w)[index]; }
  const sfloat& operator[](size_t index) const { return (&w)[index]; }

  union {
    sfloat w;
    sfloat s;
  };
  union {
    sfloat i;
    sfloat x;
  };
  union {
    sfloat j;
    sfloat y;
  };
  union {
    sfloat k;
    sfloat z;
  };
};

}  // namespace star