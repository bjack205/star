//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include <math.h>

#include "star/typedefs.h"

/*---------------------------------*/
/* Setters                         */
/*---------------------------------*/

static inline void star_SetZero4(sfloat x[4]) {
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = 0;
}

static inline void star_SetConst4(sfloat x[4], sfloat value) {
  x[0] = value;
  x[1] = value;
  x[2] = value;
  x[3] = value;
}

/*---------------------------------*/
/* Norms and Related               */
/*---------------------------------*/
static inline sfloat star_NormSquared4(const sfloat x[4]) {
  return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
}

static inline sfloat star_Norm4(const sfloat x[4]) { return sqrt(star_NormSquared4(x)); }

static inline sfloat star_OneNorm4(const sfloat x[4]) {
  return fabs(x[0]) + fabs(x[1]) + fabs(x[2]) + fabs(x[3]);
}

static inline sfloat star_InfNorm4(const sfloat x[4]) {
  sfloat x0 = fabs(x[0]);
  sfloat x1 = fabs(x[1]);
  sfloat x2 = fabs(x[2]);
  sfloat x3 = fabs(x[3]);
  sfloat norm = x1 > x0 ? x1 : x0;
  norm = x2 > norm ? x2 : norm;
  norm = x3 > norm ? x3 : norm;
  return norm;
}

static inline void star_Normalize4(sfloat out[4], const sfloat x[4]) {
  sfloat inv_norm = 1.0 / star_Norm4(x);
  out[0] *= inv_norm;
  out[1] *= inv_norm;
  out[2] *= inv_norm;
  out[3] *= inv_norm;
}

sfloat star_Dot4(const sfloat x[4], const sfloat y[4]) {
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
}

/*---------------------------------*/
/* Element-wise operations         */
/*---------------------------------*/
void star_Add4(sfloat out[4], const sfloat x[4], const sfloat y[4]) {
  out[0] = x[0] + y[0];
  out[1] = x[1] + y[1];
  out[2] = x[2] + y[2];
  out[3] = x[3] + y[3];
}

void star_Sub4(sfloat out[4], const sfloat x[4], const sfloat y[4]) {
  out[0] = x[0] - y[0];
  out[1] = x[1] - y[1];
  out[2] = x[2] - y[2];
  out[3] = x[3] - y[3];
}

void star_Mul4(sfloat out[4], const sfloat x[4], const sfloat y[4]) {
  out[0] = x[0] * y[0];
  out[1] = x[1] * y[1];
  out[2] = x[2] * y[2];
  out[3] = x[3] * y[3];
}

void star_Div4(sfloat out[4], const sfloat x[4], const sfloat y[4]) {
  out[0] = x[0] / y[0];
  out[1] = x[1] / y[1];
  out[2] = x[2] / y[2];
  out[3] = x[3] / y[3];
}

void star_UnaryMap4(sfloat out[4], const sfloat x[4], sfloat (*function)(sfloat)) {
  out[0] = function(x[0]);
  out[1] = function(x[1]);
  out[2] = function(x[2]);
  out[3] = function(x[3]);
}

void star_BinaryMap4(sfloat out[4], const sfloat x[4], const sfloat y[4],
                    sfloat (*function)(sfloat, sfloat)) {
  out[0] = function(x[0], y[0]);
  out[1] = function(x[1], y[1]);
  out[2] = function(x[2], y[2]);
  out[3] = function(x[3], y[3]);
}
