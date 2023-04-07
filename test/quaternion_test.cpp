//
// Created by Brian Jackson on 4/1/23.
// Copyright (c) 2023. All rights reserved.
//

#include <gtest/gtest.h>
#include <math.h>

extern "C" {
#include "star/matrix3.h"
#include "star/quaternion.h"
}

#define EPS 1e-8

TEST(QuaternionTest, Norm) {
  double q[4] = {1, 2, 3, 4};
  EXPECT_DOUBLE_EQ(5.4772255750516612, star_QuatNorm(q));
  EXPECT_DOUBLE_EQ(30, star_QuatNormSquared(q));
}

TEST(QuaternionTest, QuatVecNorm) {
  double q[4] = {1, 2, 3, 4};
  EXPECT_DOUBLE_EQ(5.385164807134504, star_QuatVecNorm(q));
  EXPECT_DOUBLE_EQ(29, star_QuatVecNormSquared(q));
}

TEST(QuaternionTest, Identity) {
  double q[4] = {1, 2, 3, 4};
  star_QuatIdentity(q);
  EXPECT_DOUBLE_EQ(1, q[0]);
  EXPECT_DOUBLE_EQ(0, q[1]);
  EXPECT_DOUBLE_EQ(0, q[2]);
  EXPECT_DOUBLE_EQ(0, q[3]);
  EXPECT_DOUBLE_EQ(1, star_QuatNorm(q));
}

TEST(QuaternionTest, Normalize) {
  double q[4] = {1, 2, 3, 4};
  double q_normalized[4];
  double q_norm = star_QuatNorm(q);
  star_QuatNormalize(q_normalized, q);
  EXPECT_DOUBLE_EQ(1, star_QuatNorm(q_normalized));
  EXPECT_DOUBLE_EQ(q[0] / q_norm, q_normalized[0]);
  EXPECT_DOUBLE_EQ(q[1] / q_norm, q_normalized[1]);
  EXPECT_DOUBLE_EQ(q[2] / q_norm, q_normalized[2]);
  EXPECT_DOUBLE_EQ(q[3] / q_norm, q_normalized[3]);
}

TEST(QuaternionTest, NormalizeAliased) {
  double q[4] = {1, 2, 3, 4};
  double q_norm = star_QuatNorm(q);
  star_QuatNormalize(q, q);
  EXPECT_DOUBLE_EQ(1, star_QuatNorm(q));
  EXPECT_DOUBLE_EQ(1 / q_norm, q[0]);
  EXPECT_DOUBLE_EQ(2 / q_norm, q[1]);
  EXPECT_DOUBLE_EQ(3 / q_norm, q[2]);
  EXPECT_DOUBLE_EQ(4 / q_norm, q[3]);
}

TEST(QuaternionTest, Flip) {
  double q[4] = {1, 2, 3, 4};
  double q_flip[4];
  star_QuatFlip(q_flip, q);
  EXPECT_DOUBLE_EQ(-1, q_flip[0]);
  EXPECT_DOUBLE_EQ(-2, q_flip[1]);
  EXPECT_DOUBLE_EQ(-3, q_flip[2]);
  EXPECT_DOUBLE_EQ(-4, q_flip[3]);
}

TEST(QuaternionTest, QuatVec) {
  double q[4] = {1, 2, 3, 4};
  double vec[3];
  star_QuatVec(vec, q);
  EXPECT_DOUBLE_EQ(2, vec[0]);
  EXPECT_DOUBLE_EQ(3, vec[1]);
  EXPECT_DOUBLE_EQ(4, vec[2]);
}

TEST(QuaternionTest, Conjugate) {
  double q[4] = {1, 2, 3, 4};
  double q_conj[4];
  star_QuatConjugate(q_conj, q);
  EXPECT_DOUBLE_EQ(1, q_conj[0]);
  EXPECT_DOUBLE_EQ(-2, q_conj[1]);
  EXPECT_DOUBLE_EQ(-3, q_conj[2]);
  EXPECT_DOUBLE_EQ(-4, q_conj[3]);
}

TEST(QuaternionTest, QuatInverse) {
  double q[4] = {1, 2, 3, 4};
  double q_inv[4];
  star_QuatInverse(q_inv, q);
  EXPECT_DOUBLE_EQ(1.0 / 30.0, q_inv[0]);
  EXPECT_DOUBLE_EQ(-2.0 / 30.0, q_inv[1]);
  EXPECT_DOUBLE_EQ(-3.0 / 30.0, q_inv[2]);
  EXPECT_DOUBLE_EQ(-4.0 / 30.0, q_inv[3]);
}

TEST(QuaternionTest, QuatCompose) {
  double q1[4] = {1, 2, 3, 4};
  double q2[4] = {5, 6, 7, 8};
  double q3[4];
  star_QuatCompose(q3, q1, q2);
  EXPECT_DOUBLE_EQ(-60, q3[0]);
  EXPECT_DOUBLE_EQ(12, q3[1]);
  EXPECT_DOUBLE_EQ(30, q3[2]);
  EXPECT_DOUBLE_EQ(24, q3[3]);
}

TEST(QuaternionTest, ComposeInverse) {
  double q[4] = {1, 2, 3, 4};
  double q_inv[4];
  double q_compose[4];
  star_QuatInverse(q_inv, q);
  star_QuatCompose(q_compose, q, q_inv);
  EXPECT_NEAR(1, q_compose[0], EPS);
  EXPECT_NEAR(0, q_compose[1], EPS);
  EXPECT_NEAR(0, q_compose[2], EPS);
  EXPECT_NEAR(0, q_compose[3], EPS);
}

TEST(QuaternionTest, ComposeIdentity) {
  double q[4] = {1, 2, 3, 4};
  double qI[4];
  double q_compose[4];
  star_QuatIdentity(qI);
  star_QuatCompose(q_compose, q, qI);
  EXPECT_DOUBLE_EQ(1, q_compose[0]);
  EXPECT_DOUBLE_EQ(2, q_compose[1]);
  EXPECT_DOUBLE_EQ(3, q_compose[2]);
  EXPECT_DOUBLE_EQ(4, q_compose[3]);

  star_QuatCompose(q_compose, qI, q);
  EXPECT_DOUBLE_EQ(1, q_compose[0]);
  EXPECT_DOUBLE_EQ(2, q_compose[1]);
  EXPECT_DOUBLE_EQ(3, q_compose[2]);
  EXPECT_DOUBLE_EQ(4, q_compose[3]);
}

TEST(QuaternionTest, ComposeLeft) {
  double q1[4] = {1, 2, 3, 4};
  double q2[4] = {5, 6, 7, 8};
  double q3[4];
  star_QuatComposeLeft(q3, q2, q1);
  EXPECT_DOUBLE_EQ(-60, q3[0]);
  EXPECT_DOUBLE_EQ(12, q3[1]);
  EXPECT_DOUBLE_EQ(30, q3[2]);
  EXPECT_DOUBLE_EQ(24, q3[3]);
}

TEST(QuaternionTest, Log) {
  double q[4] = {1, 2, 3, 4};
  double q_log[4];
  double phi[3] = {0.515190292664085, 0.7727854389961275, 1.03038058532817};
  star_QuatLog(q_log, q);
  EXPECT_DOUBLE_EQ(log(star_QuatNorm(q)), q_log[0]);
  EXPECT_DOUBLE_EQ(phi[0], q_log[1]);
  EXPECT_DOUBLE_EQ(phi[1], q_log[2]);
  EXPECT_DOUBLE_EQ(phi[2], q_log[3]);

  // Normalize and check that the log is zero.
  star_QuatNormalize(q, q);
  star_QuatLog(q_log, q);
  EXPECT_NEAR(0, q_log[0], 1e-15);
  EXPECT_DOUBLE_EQ(phi[0], q_log[1]);
  EXPECT_DOUBLE_EQ(phi[1], q_log[2]);
  EXPECT_DOUBLE_EQ(phi[2], q_log[3]);
}

TEST(QuaternionTest, Exp) {
  double u[3] = {0.2672612419124244, 0.5345224838248488, 0.8017837257372732};
  double theta = 3.2;
  double phi[3] = {theta * u[0], theta * u[1], theta * u[2]};
  double q_expected[4] = {cos(theta), sin(theta) * u[0], sin(theta) * u[1],
                          sin(theta) * u[2]};
  double q[4] = {0, phi[0], phi[1], phi[2]};
  double q_exp[4];
  star_QuatExp(q_exp, q);
  EXPECT_NEAR(q_expected[0], q_exp[0], EPS);
  EXPECT_NEAR(q_expected[1], q_exp[1], EPS);
  EXPECT_NEAR(q_expected[2], q_exp[2], EPS);
  EXPECT_NEAR(q_expected[3], q_exp[3], EPS);
  // Check Unit Norm
  EXPECT_NEAR(1, star_QuatNorm(q_exp), 1e-15);

  // Check with non-zero scalar
  double scalar = 1.2;
  q[0] = scalar;
  star_QuatExp(q_exp, q);
  EXPECT_NEAR(q_expected[0] * exp(scalar), q_exp[0], EPS);
  EXPECT_NEAR(q_expected[1] * exp(scalar), q_exp[1], EPS);
  EXPECT_NEAR(q_expected[2] * exp(scalar), q_exp[2], EPS);
  EXPECT_NEAR(q_expected[3] * exp(scalar), q_exp[3], EPS);
}

TEST(QuaternionTest, Exp_SmallAngle) {
  double u[3] = {0.2672612419124244, 0.5345224838248488, 0.8017837257372732};
  double theta = 1e-8;
  double phi[3] = {theta * u[0], theta * u[1], theta * u[2]};
  double q_expected[4] = {1, theta * u[0], theta * u[1], theta * u[2]};
  double q[4] = {0, phi[0], phi[1], phi[2]};
  double q_exp[4];
  star_QuatExp(q_exp, q);
  EXPECT_NEAR(q_expected[0], q_exp[0], EPS);
  EXPECT_NEAR(q_expected[1], q_exp[1], EPS);
  EXPECT_NEAR(q_expected[2], q_exp[2], EPS);
  EXPECT_NEAR(q_expected[3], q_exp[3], EPS);

  // Check Unit Norm
  EXPECT_NEAR(1, star_QuatNorm(q_exp), 1e-15);
}

TEST(QuaternionTest, LogExp) {
  double q[4] = {1, 2, 3, 4};
  double q_log[4];
  double q_exp[4];
  star_QuatLog(q_log, q);
  star_QuatExp(q_exp, q_log);
  EXPECT_NEAR(q[0], q_exp[0], EPS);
  EXPECT_NEAR(q[1], q_exp[1], EPS);
  EXPECT_NEAR(q[2], q_exp[2], EPS);
  EXPECT_NEAR(q[3], q_exp[3], EPS);
}

TEST(QuaternionTest, ExpLog) {
  double q[4] = {1, 2, 3, 4};
  double q_log[4];
  double q_exp[4];
  star_QuatNormalize(q, q);
  star_QuatExp(q_exp, q);
  star_QuatLog(q_log, q_exp);
  EXPECT_NEAR(q[0], q_log[0], EPS);
  EXPECT_NEAR(q[1], q_log[1], EPS);
  EXPECT_NEAR(q[2], q_log[2], EPS);
  EXPECT_NEAR(q[3], q_log[3], EPS);
}

TEST(QuaternionTest, RotateActive) {
  // Rotate 90 degrees about Z
  double phi[3] = {0, 0, M_PI_2};
  double q[4];
  star_QuatExpm(q, phi);
  EXPECT_NEAR(sqrt(2) / 2, q[0], EPS);
  EXPECT_NEAR(0, q[1], EPS);
  EXPECT_NEAR(0, q[2], EPS);
  EXPECT_NEAR(sqrt(2) / 2, q[3], EPS);

  double v[3] = {1, -2, 3};
  double v_rot[3];
  star_QuatRotateActive(v_rot, q, v);
  EXPECT_DOUBLE_EQ(-v[1], v_rot[0]);
  EXPECT_DOUBLE_EQ(+v[0], v_rot[1]);
  EXPECT_DOUBLE_EQ(+v[2], v_rot[2]);

  // Rotate 90 degrees about X
  phi[0] = M_PI_2;
  phi[1] = 0;
  phi[2] = 0;
  star_QuatExpm(q, phi);
  star_QuatRotateActive(v_rot, q, v);
  EXPECT_DOUBLE_EQ(+v[0], v_rot[0]);
  EXPECT_DOUBLE_EQ(-v[2], v_rot[1]);
  EXPECT_DOUBLE_EQ(+v[1], v_rot[2]);

  // Rotate 90 degrees about Y
  phi[0] = 0;
  phi[1] = M_PI_2;
  phi[2] = 0;
  star_QuatExpm(q, phi);
  star_QuatRotateActive(v_rot, q, v);
  EXPECT_NEAR(+v[2], v_rot[0], EPS);
  EXPECT_NEAR(+v[1], v_rot[1], EPS);
  EXPECT_NEAR(-v[0], v_rot[2], EPS);
}

TEST(QuaternionTest, RotatePassive) {
  // Rotate 90 degrees about Z
  double phi[3] = {0, 0, M_PI_2};
  double q[4];
  star_QuatExpm(q, phi);
  EXPECT_NEAR(sqrt(2) / 2, q[0], EPS);
  EXPECT_NEAR(0, q[1], EPS);
  EXPECT_NEAR(0, q[2], EPS);
  EXPECT_NEAR(sqrt(2) / 2, q[3], EPS);

  double v[3] = {1, -2, 3};
  double v_rot[3];
  star_QuatRotatePassive(v_rot, q, v);
  EXPECT_DOUBLE_EQ(+v[1], v_rot[0]);
  EXPECT_DOUBLE_EQ(-v[0], v_rot[1]);
  EXPECT_DOUBLE_EQ(+v[2], v_rot[2]);

  // Rotate 90 degrees about X
  phi[0] = M_PI_2;
  phi[1] = 0;
  phi[2] = 0;
  star_QuatExpm(q, phi);
  star_QuatRotatePassive(v_rot, q, v);
  EXPECT_DOUBLE_EQ(+v[0], v_rot[0]);
  EXPECT_DOUBLE_EQ(+v[2], v_rot[1]);
  EXPECT_DOUBLE_EQ(-v[1], v_rot[2]);

  // Rotate 90 degrees about Y
  phi[0] = 0;
  phi[1] = M_PI_2;
  phi[2] = 0;
  star_QuatExpm(q, phi);
  star_QuatRotatePassive(v_rot, q, v);
  EXPECT_NEAR(-v[2], v_rot[0], EPS);
  EXPECT_NEAR(+v[1], v_rot[1], EPS);
  EXPECT_NEAR(+v[0], v_rot[2], EPS);
}

TEST(QuaternionTest, QuatPure) {
  double v[3] = {1, -2, 3};
  double q_pure[4];
  star_QuatPure(q_pure, v);
  EXPECT_DOUBLE_EQ(0, q_pure[0]);
  EXPECT_DOUBLE_EQ(v[0], q_pure[1]);
  EXPECT_DOUBLE_EQ(v[1], q_pure[2]);
  EXPECT_DOUBLE_EQ(v[2], q_pure[3]);
}

TEST(QuaternionTest, QuatPureCompose) {
  double q[4] = {1, 2, 3, 4};
  double v[3] = {1, -2, 3};
  double q_pure[4];
  double q_composed[4];
  double q_composed_pure[4];
  star_QuatPure(q_pure, v);
  star_QuatCompose(q_composed, q, q_pure);
  star_QuatComposePure(q_composed_pure, q, v);
  EXPECT_NEAR(q_composed[0], q_composed_pure[0], EPS);
  EXPECT_NEAR(q_composed[1], q_composed_pure[1], EPS);
  EXPECT_NEAR(q_composed[2], q_composed_pure[2], EPS);
  EXPECT_NEAR(q_composed[3], q_composed_pure[3], EPS);
}

TEST(QuaternionConversions, ToRotMatActive) {
  double q[4] = {1, 2, 3, 4};
  star_QuatNormalize(q, q);
  double R[9];
  star_QuatToRotMatActive(R, q);
  double x[3] = {1, 2, 3};
  double y1[3];
  double y2[3];
  star_QuatRotateActive(y1, q, x);
  star_VecMul33(y2, R, x);
  EXPECT_NEAR(y1[0], y2[0], EPS);
  EXPECT_NEAR(y1[1], y2[1], EPS);
  EXPECT_NEAR(y1[2], y2[2], EPS);
}

TEST(QuaternionConversions, ToRotMatPassive) {
  double q[4] = {1, 2, 3, 4};
  star_QuatNormalize(q, q);
  double R[9];
  star_QuatToRotMatPassive(R, q);
  double x[3] = {1, 2, 3};
  double y1[3];
  double y2[3];
  star_QuatRotatePassive(y1, q, x);
  star_VecMul33(y2, R, x);
  EXPECT_NEAR(y1[0], y2[0], EPS);
  EXPECT_NEAR(y1[1], y2[1], EPS);
  EXPECT_NEAR(y1[2], y2[2], EPS);
}

TEST(QuaternionConversions, ToRodriguesParam) {
  double q[4] = {1, 2, 3, 4};
  star_QuatNormalize(q, q);
  double phi[3];

  star_QuatToRodriguesParam(phi, q);
  double q2[4];
  star_RodriguesParamToQuat(q2, phi);
  EXPECT_NEAR(q[0], q2[0], EPS);
  EXPECT_NEAR(q[1], q2[1], EPS);
  EXPECT_NEAR(q[2], q2[2], EPS);
  EXPECT_NEAR(q[3], q2[3], EPS);
}

TEST(QuaternionConversions, ToMRP) {
  double q[4] = {1, 2, 3, 4};
  star_QuatNormalize(q, q);
  double mrp[3];

  star_QuatToMRP(mrp, q);
  double q2[4];
  star_MRPToQuat(q2, mrp);
  EXPECT_NEAR(q[0], q2[0], EPS);
  EXPECT_NEAR(q[1], q2[1], EPS);
  EXPECT_NEAR(q[2], q2[2], EPS);
  EXPECT_NEAR(q[3], q2[3], EPS);
}

TEST(QuaternionConversions, ToEuler123) {
  double q[4] = {1, 2, 3, 4};
  star_QuatNormalize(q, q);
  double euler[3];
  double euler_expected[3] = {-1.373400766945016, 0.8232119771258758, -2.9441970937399127};

  star_QuatToEulerXYZ(euler, q);
  EXPECT_NEAR(euler_expected[0], euler[0], EPS);
  EXPECT_NEAR(euler_expected[1], euler[1], EPS);
  EXPECT_NEAR(euler_expected[2], euler[2], EPS);

  // Rotate about x
  double angle = 0.2;
  double phi[3] = {angle, 0, 0};
  star_QuatExpm(q, phi);
  star_QuatToEulerXYZ(euler, q);
  EXPECT_NEAR(euler[0], angle, EPS);
  EXPECT_NEAR(euler[1], 0, EPS);
  EXPECT_NEAR(euler[2], 0, EPS);

  // Rotate about y
  phi[0] = 0;
  phi[1] = angle;
  star_QuatExpm(q, phi);
  star_QuatToEulerXYZ(euler, q);
  EXPECT_NEAR(euler[0], 0, EPS);
  EXPECT_NEAR(euler[1], angle, EPS);
  EXPECT_NEAR(euler[2], 0, EPS);

  // Rotate about z
  phi[1] = 0;
  phi[2] = angle;
  star_QuatExpm(q, phi);
  star_QuatToEulerXYZ(euler, q);
  EXPECT_NEAR(euler[0], 0, EPS);
  EXPECT_NEAR(euler[1], 0, EPS);
  EXPECT_NEAR(euler[2], angle, EPS);
}

TEST(QuaternionConversions, ToEuler321) {
  double q[4] = {1, 2, 3, 4};
  star_QuatNormalize(q, q);
  double euler[3];
  double euler_expected[3] = {2.356194490192345, -0.3398369094541219, 1.4288992721907328};

  star_QuatToEulerZYX(euler, q);
  EXPECT_NEAR(euler_expected[0], euler[0], EPS);
  EXPECT_NEAR(euler_expected[1], euler[1], EPS);
  EXPECT_NEAR(euler_expected[2], euler[2], EPS);

  // Rotate about x
  double angle = 0.2;
  double phi[3] = {angle, 0, 0};
  star_QuatExpm(q, phi);
  star_QuatToEulerZYX(euler, q);
  EXPECT_NEAR(euler[0], 0, EPS);
  EXPECT_NEAR(euler[1], 0, EPS);
  EXPECT_NEAR(euler[2], angle, EPS);

  // Rotate about y
  phi[0] = 0;
  phi[1] = angle;
  star_QuatExpm(q, phi);
  star_QuatToEulerZYX(euler, q);
  EXPECT_NEAR(euler[0], 0, EPS);
  EXPECT_NEAR(euler[1], angle, EPS);
  EXPECT_NEAR(euler[2], 0, EPS);

  // Rotate about z
  phi[1] = 0;
  phi[2] = angle;
  star_QuatExpm(q, phi);
  star_QuatToEulerZYX(euler, q);
  EXPECT_NEAR(euler[0], angle, EPS);
  EXPECT_NEAR(euler[1], 0, EPS);
  EXPECT_NEAR(euler[2], 0, EPS);
}