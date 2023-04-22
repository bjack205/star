//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#include <fmt/core.h>
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "star/Quaternion.hpp"
#include "star/Vec3.hpp"

#define EPS 1e-8

using namespace star;

TEST(QuaternionClass, Norm) {
  Quaternion q(1, 2, 3, 4);
  EXPECT_FLOAT_EQ(q.Norm(), sqrt(1 + 4 + 9 + 16));
}

TEST(QuaternionClass, VecNorm) {
  Quaternion q(1, 2, 3, 4);
  EXPECT_FLOAT_EQ(q.VecNorm(), sqrt(4 + 9 + 16));
}

TEST(QuaternionClass, VecNormSquared) {
  Quaternion q(1, 2, 3, 4);
  EXPECT_FLOAT_EQ(q.VecNormSquared(), 4 + 9 + 16);
}

TEST(QuaternionClass, NormalizeInPlace) {
  Quaternion q(1, 2, 3, 4);
  q.NormalizeInPlace();
  EXPECT_FLOAT_EQ(q.Norm(), 1);
}

TEST(QuaternionClass, ImplicitConversion) {
  Quaternion q(1, 2, 3, 4);
  Vec4 v = q;
  EXPECT_FLOAT_EQ(v[0], q[0]);
  EXPECT_FLOAT_EQ(v[1], q[1]);
  EXPECT_FLOAT_EQ(v[2], q[2]);
  EXPECT_FLOAT_EQ(v[3], q[3]);
}

TEST(QuaternionClass, FromVector) {
  std::vector<sfloat> v = {1, 2, 3, 4};
  Quaternion q(v);
  EXPECT_FLOAT_EQ(q[0], v[0]);
  EXPECT_FLOAT_EQ(q[1], v[1]);
  EXPECT_FLOAT_EQ(q[2], v[2]);
  EXPECT_FLOAT_EQ(q[3], v[3]);
}

TEST(QuaternionClass, Normalize) {
  Quaternion q(1, 2, 3, 4);
  Quaternion q2 = q.Normalize();
  EXPECT_FLOAT_EQ(q2.Norm(), 1);
}

TEST(QuaternionClass, Indentity) {
  Quaternion q = Quaternion::Identity();
  EXPECT_FLOAT_EQ(q.Norm(), 1);
  EXPECT_FLOAT_EQ(q[0], 1);
  EXPECT_FLOAT_EQ(q[1], 0);
  EXPECT_FLOAT_EQ(q[2], 0);
  EXPECT_FLOAT_EQ(q[3], 0);
}

TEST(QuaternionClass, Pure) {
  Quaternion q = Quaternion::Pure(1, 2, 3);
  EXPECT_FLOAT_EQ(q.Norm(), sqrt(1 + 4 + 9));
  EXPECT_FLOAT_EQ(q[0], 0);
  EXPECT_FLOAT_EQ(q[1], 1);
  EXPECT_FLOAT_EQ(q[2], 2);
  EXPECT_FLOAT_EQ(q[3], 3);

  // From Vec3
  Vec3 v(1, 2, 3);
  Quaternion q2 = Quaternion::Pure(v);
  EXPECT_FLOAT_EQ(q2.Norm(), sqrt(1 + 4 + 9));
  EXPECT_FLOAT_EQ(q2[0], 0);
  EXPECT_FLOAT_EQ(q2[1], 1);
  EXPECT_FLOAT_EQ(q2[2], 2);
  EXPECT_FLOAT_EQ(q2[3], 3);
}

TEST(QuaternionClass, Expm) {
  // Rotate 90 degrees about the x-axis
  Quaternion q = Quaternion::Expm(M_PI / 2, 0, 0);
  EXPECT_FLOAT_EQ(q.Norm(), 1);
  EXPECT_NEAR(q[0], sqrt(2) / 2.0, EPS);
  EXPECT_NEAR(q[1], sqrt(2) / 2.0, EPS);
  EXPECT_NEAR(q[2], 0, EPS);
  EXPECT_NEAR(q[3], 0, EPS);

  Quaternion q2 = Quaternion::RotX(M_PI / 2);
  EXPECT_TRUE(q.IsApprox(q2));

  // Rotate 45 degrees about the y-axis
  sfloat angle = M_PI / 4;
  q = Quaternion::Expm(0, angle, 0);
  EXPECT_FLOAT_EQ(q2.Norm(), 1);
  EXPECT_NEAR(q[0], cos(angle / 2), EPS);
  EXPECT_NEAR(q[1], 0, EPS);
  EXPECT_NEAR(q[2], sin(angle / 2), EPS);
  EXPECT_NEAR(q[3], 0, EPS);

  q2 = Quaternion::RotY(angle);
  EXPECT_TRUE(q.IsApprox(q2));

  // Rotate 30 degrees about the z-axis
  angle = M_PI / 6;
  q = Quaternion::Expm(0, 0, angle);
  EXPECT_FLOAT_EQ(q.Norm(), 1);
  EXPECT_NEAR(q[0], cos(angle / 2), EPS);
  EXPECT_NEAR(q[1], 0, EPS);
  EXPECT_NEAR(q[2], 0, EPS);
  EXPECT_NEAR(q[3], sin(angle / 2), EPS);

  q2 = Quaternion::RotZ(angle);
  EXPECT_TRUE(q.IsApprox(q2));
}

TEST(QuaternionClass, AngleBetween) {
  Quaternion q = {1, 2, 3, 4};
  q.NormalizeInPlace();
  sfloat angle = 1.2;
  Vec3 phi = Vec3(1, 2, 3).Normalize();
  phi *= angle;
  Quaternion q2 = q.Compose(Quaternion::Expm(phi));
  EXPECT_FLOAT_EQ(q.AngleBetween(q2), angle);
}
