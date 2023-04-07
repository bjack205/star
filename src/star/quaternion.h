//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "typedefs.h"

// Scalar values
double star_QuatNorm(const double q[4]);
double star_QuatNormSquared(const double q[4]);
double star_QuatVecNorm(const double q[4]);
double star_QuatVecNormSquared(const double q[4]);

// Quaternion operations
void star_QuatIdentity(double q[4]);
void star_QuatNormalize(double q_normalized[4], double q[4]);
void star_QuatFlip(double q_flip[4], const double q[4]);
void star_QuatVec(double vec[3], const double q[4]);
void star_QuatConjugate(double q_conj[4], const double q[4]);
void star_QuatInverse(double q_inv[4], const double q[4]);
void star_QuatCompose(double q12[4], const double q1[4], const double q2[4]);
void star_QuatComposeLeft(double q21[4], const double q1[4], const double q2[4]);

// Operations on vectors
void star_QuatLogm(double phi[4], const double q[4]);
void star_QuatLog(double q_log[4], const double q[4]);
void star_QuatExpm(double q[4], const double phi[3]);
void star_QuatExp(double q_exp[4], const double q[3]);
void star_QuatRotateActive(double v_rot[4], const double q[4], const double v[3]);
void star_QuatRotatePassive(double* v_rot, const double q[3], const double v[3]);
void star_QuatPure(double q[4], const double x[3]);
void star_QuatComposePure(double qv[4], const double q[4], const double v[3]);

// Conversions
void star_QuatToRotMatActive(double Q[9], const double q[4]);
void star_QuatToRotMatPassive(double Q[9], const double q[4]);

void star_QuatToRodriguesParam(double g[3], const double q[4]);
void star_RodriguesParamToQuat(double q[4], const double g[3]);

void star_QuatToMRP(double p[3], const double q[4]);
void star_MRPToQuat(double q[4], const double p[3]);

void star_QuatToAxisAngle(double aa[4], const double q[4]);
void star_QuatToEulerXYZ(double e[3], const double q[4]);
void star_QuatToEulerZYX(double e[3], const double q[4]);

// Jacobians
void star_QuatRotateActiveJacobian(double* D, const double* q, const double* x);
void star_QuatRotatePassiveJacobian(double* D, const double* q, const double* x);

// Matrices
void star_SkewSymmetricMatrix(double* S, const double* x);
void star_LMat(double* L, const double* q);
void star_RMat(double* R, const double* q);
void star_GMat(double* G, double* q);



void qmat_err(double* phi, const double* q1, const double* q2);
void qmat_adderr(double* q, const double* q0, const double* phi);

void qmat_cay(double* q, const double* phi);
void qmat_icay(double* phi, const double* q);
void qmat_dcay(double* D, const double* phi);
