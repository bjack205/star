//
// Created by Brian Jackson on 4/25/23.
// Copyright (c) 2023. All rights reserved.
//

#include <gtest/gtest.h>

#include "star/RotMat.hpp"

using namespace star;

TEST(RotMat, DefaultConstructor) {
  RotMat<Active> R;
  EXPECT_FLOAT_EQ(R[0], 1);
  EXPECT_FLOAT_EQ(R[1], 0);
  EXPECT_FLOAT_EQ(R[2], 0);
  EXPECT_FLOAT_EQ(R[3], 0);
  EXPECT_FLOAT_EQ(R[4], 1);
  EXPECT_FLOAT_EQ(R[5], 0);
  EXPECT_FLOAT_EQ(R[6], 0);
  EXPECT_FLOAT_EQ(R[7], 0);
  EXPECT_FLOAT_EQ(R[8], 1);
}

TEST(RotMat, RotX) {
  sfloat angle = 0.5;
  RotMat<Active> R = RotMat<Active>::RotX(angle);
  EXPECT_FLOAT_EQ(R(0, 0), 1);
  EXPECT_FLOAT_EQ(R(0, 1), 0);
  EXPECT_FLOAT_EQ(R(0, 2), 0);
  EXPECT_FLOAT_EQ(R(1, 0), 0);
  EXPECT_FLOAT_EQ(R(1, 1), cos(angle));
  EXPECT_FLOAT_EQ(R(1, 2), -sin(angle));
  EXPECT_FLOAT_EQ(R(2, 0), 0);
  EXPECT_FLOAT_EQ(R(2, 1), sin(angle));
  EXPECT_FLOAT_EQ(R(2, 2), cos(angle));

  RotMat<Passive> A = RotMat<Passive>::RotX(angle);
  EXPECT_FLOAT_EQ(A(0, 0), 1);
  EXPECT_FLOAT_EQ(A(0, 1), 0);
  EXPECT_FLOAT_EQ(A(0, 2), 0);
  EXPECT_FLOAT_EQ(A(1, 0), 0);
  EXPECT_FLOAT_EQ(A(1, 1), cos(angle));
  EXPECT_FLOAT_EQ(A(1, 2), sin(angle));
  EXPECT_FLOAT_EQ(A(2, 0), 0);
  EXPECT_FLOAT_EQ(A(2, 1), -sin(angle));
  EXPECT_FLOAT_EQ(A(2, 2), cos(angle));
}

TEST(RotMat, Transpose) {
  RotMat<Active> R = RotMat<Active>::RotX(0.5);
  RotMat<Active> R_T = R.Transpose();
  RotMat<Passive> A = RotMat<Passive>::RotX(0.5);
  for (int i = 0; i < 9; i++) {
    EXPECT_FLOAT_EQ(R_T[i], A[i]);
  }
}