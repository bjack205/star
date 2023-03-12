//
// Created by Brian Jackson on 3/9/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <math.h>
#include <stdint.h>

#include "typedefs.h"

/*---------------------------------*/
/* Setters                         */
/*---------------------------------*/

enum star_Axis { star_Xaxis, star_Yaxis, star_Zaxis };

static inline void star_SetZero(sfloat x[3]) {
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
}

static inline void star_SetConst(sfloat x[3], sfloat value) {
  x[0] = value;
  x[1] = value;
  x[2] = value;
}

static inline void star_SetAxis(sfloat x[3], enum star_Axis ax) {
  switch (ax) {
    case star_Xaxis:
      x[0] = 1;
      x[1] = 0;
      x[2] = 0;
      break;
    case star_Yaxis:
      x[0] = 0;
      x[1] = 1;
      x[2] = 0;
      break;
    case star_Zaxis:
      x[0] = 0;
      x[1] = 0;
      x[2] = 1;
      break;
    default:
      x[0] = 0;
      x[1] = 0;
      x[2] = 0;
  }
}

/*---------------------------------*/
/* Norms and Related               */
/*---------------------------------*/
static inline sfloat star_NormSquared3(const sfloat x[3]) {
  return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

static inline sfloat star_Norm3(const sfloat x[3]) { return sqrt(star_NormSquared3(x)); }

static inline sfloat star_OneNorm3(const sfloat x[3]) {
  return fabs(x[0]) + fabs(x[1]) + fabs(x[2]);
}

static inline sfloat star_InfNorm3(const sfloat x[3]) {
  sfloat x0 = fabs(x[0]);
  sfloat x1 = fabs(x[1]);
  sfloat x2 = fabs(x[2]);
  sfloat norm = x1 > x0 ? x1 : x0;
  norm = x2 > norm ? x2 : norm;
  return norm;
}

static inline void star_Normalize3(sfloat out[3], const sfloat x[3]) {
  sfloat inv_norm = 1.0 / star_Norm3(x);
  out[0] *= inv_norm;
  out[1] *= inv_norm;
  out[2] *= inv_norm;
}

sfloat star_Dot3(const sfloat x[3], const sfloat y[3]) {
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

/*---------------------------------*/
/* Element-wise operations         */
/*---------------------------------*/
void star_Add3(sfloat out[3], const sfloat x[3], const sfloat y[3]) {
  out[0] = x[0] + y[0];
  out[1] = x[1] + y[1];
  out[2] = x[2] + y[2];
}

void star_Sub3(sfloat out[3], const sfloat x[3], const sfloat y[3]) {
  out[0] = x[0] - y[0];
  out[1] = x[1] - y[1];
  out[2] = x[2] - y[2];
}

void star_Mul3(sfloat out[3], const sfloat x[3], const sfloat y[3]) {
  out[0] = x[0] * y[0];
  out[1] = x[1] * y[1];
  out[2] = x[2] * y[2];
}

void star_Div3(sfloat out[3], const sfloat x[3], const sfloat y[3]) {
  out[0] = x[0] / y[0];
  out[1] = x[1] / y[1];
  out[2] = x[2] / y[2];
}

void star_UnaryMap(sfloat out[3], const sfloat x[3], sfloat (*function)(sfloat)) {
  out[0] = function(x[0]);
  out[1] = function(x[1]);
  out[2] = function(x[2]);
}

void star_BinaryMap(sfloat out[3], const sfloat x[3], const sfloat y[3],
                    sfloat (*function)(sfloat, sfloat)) {
  out[0] = function(x[0], y[0]);
  out[1] = function(x[1], y[1]);
  out[2] = function(x[2], y[2]);
}
