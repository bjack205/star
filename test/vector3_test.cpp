//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#include <cmath>
#include <vector>
#include <limits>

#include <gtest/gtest.h>

#include "star/star.hpp"

using namespace star;

constexpr sfloat EPS = 10 * std::numeric_limits<sfloat>::epsilon();

TEST(Vector3, Constructor) {
  Vec3 x(1, 2, 3);
  EXPECT_DOUBLE_EQ(x[0], 1);
  EXPECT_DOUBLE_EQ(x[1], 2);
  EXPECT_DOUBLE_EQ(x[2], 3);

  Vec3 y{3, 4, 5};
  EXPECT_DOUBLE_EQ(y[0], 3);
  EXPECT_DOUBLE_EQ(y[1], 4);
  EXPECT_DOUBLE_EQ(y[2], 5);

  Vec3 z = {10, 11, 12};
  EXPECT_DOUBLE_EQ(z[0], 10);
  EXPECT_DOUBLE_EQ(z[1], 11);
  EXPECT_DOUBLE_EQ(z[2], 12);
}

TEST(Vector3, ConstConstructor) {
  const Vec3 x(1, 2, 3);
  EXPECT_DOUBLE_EQ(x[0], 1);
  EXPECT_DOUBLE_EQ(x[1], 2);
  EXPECT_DOUBLE_EQ(x[2], 3);

  const Vec3 y{3, 4, 5};
  EXPECT_DOUBLE_EQ(y[0], 3);
  EXPECT_DOUBLE_EQ(y[1], 4);
  EXPECT_DOUBLE_EQ(y[2], 5);

  const Vec3 z = {10, 11, 12};
  EXPECT_DOUBLE_EQ(z[0], 10);
  EXPECT_DOUBLE_EQ(z[1], 11);
  EXPECT_DOUBLE_EQ(z[2], 12);
}

TEST(Vector3, ArrayConstructor) {
  sfloat data[3] = {1, 2, 3};
  Vec3 x(data);
  EXPECT_DOUBLE_EQ(x[0], 1);
  EXPECT_DOUBLE_EQ(x[1], 2);
  EXPECT_DOUBLE_EQ(x[2], 3);
}

TEST(Vector3, VectorConstructor) {
  std::vector<sfloat> data = {1, 2, 3};
  Vec3 x(data);
  EXPECT_DOUBLE_EQ(x[0], 1);
  EXPECT_DOUBLE_EQ(x[1], 2);
  EXPECT_DOUBLE_EQ(x[2], 3);
}

TEST(Vector3, Zero) {
  Vec3 x = Vec3::Zero();
  EXPECT_DOUBLE_EQ(x[0], 0);
  EXPECT_DOUBLE_EQ(x[1], 0);
  EXPECT_DOUBLE_EQ(x[2], 0);
}

TEST(Vector3, Const) {
  const sfloat val = 1.2;
  Vec3 x = Vec3::Const(val);
  EXPECT_DOUBLE_EQ(x[0], val);
  EXPECT_DOUBLE_EQ(x[1], val);
  EXPECT_DOUBLE_EQ(x[2], val);
}

TEST(Vector3, NormSquared) {
  const Vec3 x{1, 2, 3};
  EXPECT_NEAR(x.NormSquared(), 1 + 4 + 9, EPS);
}

TEST(Vector3, Norm) {
  const Vec3 x{1, 2, 3};
  EXPECT_NEAR(x.Norm(), std::sqrt(1 + 4 + 9), EPS);
}

TEST(Vector3, OneNorm) {
  const Vec3 x{1, 2, -3};
  EXPECT_NEAR(x.OneNorm(), 1 + 2 + 3, EPS);
}

TEST(Vector3, InfNorm) {
  Vec3 x{1, 2, -3};
  EXPECT_NEAR(x.InfNorm(), 3, EPS);

  x.x = -4;  // NOLINT
  EXPECT_NEAR(x.InfNorm(), 4, EPS);

  x.y = 10;  // NOLINT
  EXPECT_NEAR(x.InfNorm(), 10, EPS);
}

TEST(Vector3, Normalize) {
  const Vec3 x = {3, 4, 5};
  sfloat norm = std::sqrt(9 + 16 + 25);   // NOLINT
  Vec3 x_bar = x.Normalize();
  EXPECT_NEAR(x.NormSquared(), 50, EPS);
  EXPECT_NEAR(x_bar[0], 3 / norm, EPS);
  EXPECT_NEAR(x_bar[1], 4 / norm, EPS);
  EXPECT_NEAR(x_bar[2], 5 / norm, EPS);
  EXPECT_NEAR(x_bar.Norm(), 1, EPS);
}

TEST(Vector3, NormalizeInPlace) {
  Vec3 x = {3, 4, 5};
  sfloat norm = std::sqrt(9 + 16 + 25);   // NOLINT
  x.NormalizeInPlace();
  EXPECT_NEAR(x[0], 3 / norm, EPS);
  EXPECT_NEAR(x[1], 4 / norm, EPS);
  EXPECT_NEAR(x[2], 5 / norm, EPS);
  EXPECT_NEAR(x.Norm(), 1, EPS);
}

TEST(Vector3, DotProduct) {
  const Vec3 x = {3, 4, 5};
  const Vec3 y = {1, 2, 3};
  sfloat dot = x.Dot(y);
  EXPECT_NEAR(dot, 3 + 8 + 15, EPS);
}

TEST(Vector3, Add) {
  const Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  Vec3 z = x.Add(y);
  EXPECT_NEAR(z[0], 2, EPS);
  EXPECT_NEAR(z[1], 6, EPS);
  EXPECT_NEAR(z[2], 8, EPS);

  Vec3 z2 = x + y;
  EXPECT_NEAR(z2[0], 2, EPS);
  EXPECT_NEAR(z2[1], 6, EPS);
  EXPECT_NEAR(z2[2], 8, EPS);
}

TEST(Vector3, Sub) {
  const Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  Vec3 z = x.Sub(y);
  EXPECT_NEAR(z[0], 4, EPS);
  EXPECT_NEAR(z[1], 2, EPS);
  EXPECT_NEAR(z[2], 2, EPS);

  Vec3 z2 = x - y;
  EXPECT_NEAR(z2[0], 4, EPS);
  EXPECT_NEAR(z2[1], 2, EPS);
  EXPECT_NEAR(z2[2], 2, EPS);
}

TEST(Vector3, Mul) {
  const Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  Vec3 z = x.Mul(y);
  EXPECT_NEAR(z[0], -3, EPS);
  EXPECT_NEAR(z[1], 8, EPS);
  EXPECT_NEAR(z[2], 15, EPS);

  Vec3 z2 = x * y;
  EXPECT_NEAR(z2[0], -3, EPS);
  EXPECT_NEAR(z2[1], 8, EPS);
  EXPECT_NEAR(z2[2], 15, EPS);
}

TEST(Vector3, Div) {
  const Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  Vec3 z = x.Div(y);
  EXPECT_NEAR(z[0], -3, EPS);
  EXPECT_NEAR(z[1], 2, EPS);
  EXPECT_NEAR(z[2], 5.0 / 3.0, EPS);

  Vec3 z2 = x / y;
  EXPECT_NEAR(z2[0], -3, EPS);
  EXPECT_NEAR(z2[1], 2, EPS);
  EXPECT_NEAR(z2[2], 5.0 / 3.0, EPS);
}

TEST(Vector3, AddInPlace) {
  Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  x.AddInPlace(y);
  EXPECT_NEAR(x[0], 2, EPS);
  EXPECT_NEAR(x[1], 6, EPS);
  EXPECT_NEAR(x[2], 8, EPS);

  x += y;
  EXPECT_NEAR(x[0], 1, EPS);
  EXPECT_NEAR(x[1], 8, EPS);
  EXPECT_NEAR(x[2], 11, EPS);
}

TEST(Vector3, SubInPlace) {
  Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  x.SubInPlace(y);
  EXPECT_NEAR(x[0], 4, EPS);
  EXPECT_NEAR(x[1], 2, EPS);
  EXPECT_NEAR(x[2], 2, EPS);

  x -= y;
  EXPECT_NEAR(x[0], 5, EPS);
  EXPECT_NEAR(x[1], 0, EPS);
  EXPECT_NEAR(x[2], -1, EPS);
}

TEST(Vector3, MulInPlace) {
  Vec3 x = {3, 4, 5};
  const Vec3 y = {-1, 2, 3};
  x.MulInPlace(y);
  EXPECT_NEAR(x[0], -3, EPS);
  EXPECT_NEAR(x[1], 8, EPS);
  EXPECT_NEAR(x[2], 15, EPS);

  x *= y;
  EXPECT_NEAR(x[0], 3, EPS);
  EXPECT_NEAR(x[1], 16, EPS);
  EXPECT_NEAR(x[2], 45, EPS);
}

TEST(Vector3, DivInPlace) {
  Vec3 x = {3, 4, 27};
  const Vec3 y = {-1, 2, 3};
  x.DivInPlace(y);
  EXPECT_NEAR(x[0], -3, EPS);
  EXPECT_NEAR(x[1], 2, EPS);
  EXPECT_NEAR(x[2], 9, EPS);

  x /= y;
  EXPECT_NEAR(x[0], 3, EPS);
  EXPECT_NEAR(x[1], 1, EPS);
  EXPECT_NEAR(x[2], 3, EPS);
}

TEST(Vector3, UnaryMap) {
  const Vec3 x = {4,9,16};
  Vec3 y = x.UnaryMap(std::sqrt);
  EXPECT_NEAR(y[0], 2, EPS);
  EXPECT_NEAR(y[1], 3, EPS);
  EXPECT_NEAR(y[2], 4, EPS);
}

sfloat foo(sfloat x, sfloat y) {
  return 2 * x - y;
}

TEST(Vector3, BinaryMap) {
  const Vec3 x = {4,9, 3};
  const Vec3 y = {-1, 3, 6};
  Vec3 z = x.BinaryMap(y, foo);
  EXPECT_NEAR(z[0], 9, EPS);
  EXPECT_NEAR(z[1], 15, EPS);
  EXPECT_NEAR(z[2], 0, EPS);
}

class MRP : public Vec3 {
 public:
  using Vec3::Vec3;
};

TEST(MRP, Constructor) {
  MRP g{1, 2, 3};
  g.Norm();
}