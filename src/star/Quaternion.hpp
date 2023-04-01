//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <cstddef>
#include <cstdint>

#include "star/Vec4.hpp"
#include "star/typedefs.h"

namespace star {

class Quaternion : public Vec4 {
 public:
  using Vec4::Vec4;

  /*---------------------------------*/
  /* Static Methods                  */
  /*---------------------------------*/
  static Quaternion Identity() { return {1,0,0,0}; };
  static Quaternion Expm(sfloat x, sfloat y, sfloat z);

  /*---------------------------------*/
  /* Mathematical operators          */
  /*---------------------------------*/
  Quaternion Exp() const;
  Quaternion Log() const;
  Quaternion Compose(const Quaternion rhs) const;
  Quaternion ComposeLeft(const Quaternion lhs) const;
};

}  // namespace star
