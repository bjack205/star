//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#include "quaternion.h"

#include <stdio.h>

#include "math.h"

double star_QuatNorm(const double* q) {
  return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}

double star_QuatNormSquared(const double* q) {
  return q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
}

double star_QuatVecNorm(const double* q) {
  return sqrt(q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}
double star_QuatVecNormSquared(const double* q) {
  return q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
}

void star_QuatIdentity(double* q) {
  q[0] = 1;
  q[1] = 0;
  q[2] = 0;
  q[3] = 0;
}

void star_QuatNormalize(double* q_normalized, double* q) {
  double n = 1 / star_QuatNorm(q);
  q_normalized[0] = q[0] * n;
  q_normalized[1] = q[1] * n;
  q_normalized[2] = q[2] * n;
  q_normalized[3] = q[3] * n;
}

void star_QuatFlip(double* q_flip, const double* q) {
  q_flip[0] = -q[0];
  q_flip[1] = -q[1];
  q_flip[2] = -q[2];
  q_flip[3] = -q[3];
}

void star_QuatVec(double* vec, const double* q) {
  vec[0] = q[1];
  vec[1] = q[2];
  vec[2] = q[3];
}

void star_QuatConjugate(double* q_conj, const double* q) {
  q_conj[0] = q[0];
  q_conj[1] = -q[1];
  q_conj[2] = -q[2];
  q_conj[3] = -q[3];
}

void star_QuatInverse(double* qinv, const double* q) {
  double n = 1 / star_QuatNormSquared(q);
  qinv[0] = q[0] * n;
  qinv[1] = -q[1] * n;
  qinv[2] = -q[2] * n;
  qinv[3] = -q[3] * n;
}

void star_QuatCompose(double* q3, const double* q1, const double* q2) {
  q3[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  q3[1] = q1[1] * q2[0] + q1[0] * q2[1] + q1[2] * q2[3] - q1[3] * q2[2];
  q3[2] = q1[2] * q2[0] + q1[3] * q2[1] + q1[0] * q2[2] - q1[1] * q2[3];
  q3[3] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];
}

void star_QuatComposeLeft(double* q3, const double* q1, const double* q2) {
  q3[0] = q2[0] * q1[0] - q2[1] * q1[1] - q2[2] * q1[2] - q2[3] * q1[3];
  q3[1] = q2[1] * q1[0] + q2[0] * q1[1] + q2[2] * q1[3] - q2[3] * q1[2];
  q3[2] = q2[2] * q1[0] + q2[3] * q1[1] + q2[0] * q1[2] - q2[1] * q1[3];
  q3[3] = q2[3] * q1[0] + q2[0] * q1[3] + q2[1] * q1[2] - q2[2] * q1[1];
}

void star_QuatLogm(double* phi, const double* q) {
  double s = q[0];
  double theta = star_QuatVecNorm(q);
  double M;
  if (theta < 1e-6) {
    if (fabs(s) < STAR_EPS) {
      phi[0] = 0;
      phi[1] = 0;
      phi[2] = 0;
      return;
    }
    M = (1 - theta * theta / (3 * s * s)) / s;
  } else {
    M = atan2(theta, s) / theta;
  }
  phi[0] = q[1] * M * 2;
  phi[1] = q[2] * M * 2;
  phi[2] = q[3] * M * 2;
}

void star_QuatLog(double* q_log, const double* q) {
  star_QuatLogm(q_log + 1, q);
  q_log[0] = log(star_QuatNorm(q));
  q_log[1] *= 0.5;
  q_log[2] *= 0.5;
  q_log[3] *= 0.5;
}

void star_QuatExpm(double* q, const double* phi) {
  double theta = sqrt(phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
  double s_theta;
  double c_theta;
  if (theta < sqrt(STAR_EPS)) {
    theta /= 2;
    s_theta = (1 - theta * theta / 6.0);
    c_theta = (1 - theta * theta / 2.0);
  } else {
    s_theta = sin(theta / 2) / theta;
    c_theta = cos(theta / 2);
  }
  q[0] = c_theta;
  q[1] = phi[0] * s_theta;
  q[2] = phi[1] * s_theta;
  q[3] = phi[2] * s_theta;
}

void star_QuatExp(double* q_exp, const double* q) {
  double phi[3] = {2 * q[1], 2 * q[2], 2 * q[3]};
  star_QuatExpm(q_exp, phi);
  double s = exp(q[0]);
  q_exp[0] *= s;
  q_exp[1] *= s;
  q_exp[2] *= s;
  q_exp[3] *= s;
}

void star_QuatRotateActive(double* v_rot, const double* q, const double* v) {
  double w = q[0];
  double x = q[1];
  double y = q[2];
  double z = q[3];
  double ww = w * w;
  double xx = x * x;
  double yy = y * y;
  double zz = z * z;
  double xy = x * y;
  double zw = z * w;
  double xz = x * z;
  double yw = y * w;
  double yz = y * z;
  double xw = x * w;

  v_rot[0] = (ww + xx - yy - zz) * v[0];
  v_rot[1] = 2 * (xy + zw) * v[0];
  v_rot[2] = 2 * (xz - yw) * v[0];
  v_rot[0] += 2 * (xy - zw) * v[1];
  v_rot[1] += (ww - xx + yy - zz) * v[1];
  v_rot[2] += 2 * (yz + xw) * v[1];
  v_rot[0] += 2 * (xz + yw) * v[2];
  v_rot[1] += 2 * (yz - xw) * v[2];
  v_rot[2] += (ww - xx - yy + zz) * v[2];
}

void star_QuatRotatePassive(double* v_rot, const double* q, const double* v) {
  double w = q[0];
  double x = q[1];
  double y = q[2];
  double z = q[3];
  double ww = w * w;
  double xx = x * x;
  double yy = y * y;
  double zz = z * z;
  double xy = x * y;
  double zw = z * w;
  double xz = x * z;
  double yw = y * w;
  double yz = y * z;
  double xw = x * w;

  v_rot[0] = (ww + xx - yy - zz) * v[0];
  v_rot[0] += 2 * (xy + zw) * v[1];
  v_rot[0] += 2 * (xz - yw) * v[2];
  v_rot[1] = 2 * (xy - zw) * v[0];
  v_rot[1] += (ww - xx + yy - zz) * v[1];
  v_rot[1] += 2 * (yz + xw) * v[2];
  v_rot[2] = 2 * (xz + yw) * v[0];
  v_rot[2] += 2 * (yz - xw) * v[1];
  v_rot[2] += (ww - xx - yy + zz) * v[2];
}

void star_QuatPure(double* q, const double* x) {
  q[0] = 0;
  q[1] = x[0];
  q[2] = x[1];
  q[3] = x[2];
}

void star_QuatComposePure(double* qv, const double* q1, const double* v) {
  qv[0] = -q1[1] * v[0] - q1[2] * v[1] - q1[3] * v[2];
  qv[1] = +q1[0] * v[0] + q1[2] * v[2] - q1[3] * v[1];
  qv[2] = +q1[3] * v[0] + q1[0] * v[1] - q1[1] * v[2];
  qv[3] = +q1[0] * v[2] + q1[1] * v[1] - q1[2] * v[0];
}

/////////////////////////////////////////////
// Conversions
/////////////////////////////////////////////

void star_QuatToRotMatActive(double* Q, const double* q) {
  double w = q[0];
  double x = q[1];
  double y = q[2];
  double z = q[3];
  double ww = w * w;
  double xx = x * x;
  double yy = y * y;
  double zz = z * z;
  double xy = x * y;
  double zw = z * w;
  double xz = x * z;
  double yw = y * w;
  double yz = y * z;
  double xw = x * w;

  Q[0 + 0] = ww + xx - yy - zz;
  Q[0 + 1] = 2 * (xy + zw);
  Q[0 + 2] = 2 * (xz - yw);
  Q[3 + 0] = 2 * (xy - zw);
  Q[3 + 1] = ww - xx + yy - zz;
  Q[3 + 2] = 2 * (yz + xw);
  Q[6 + 0] = 2 * (xz + yw);
  Q[6 + 1] = 2 * (yz - xw);
  Q[6 + 2] = ww - xx - yy + zz;
}

void star_QuatToRotMatPassive(double* Q, const double* q) {
  double w = q[0];
  double x = q[1];
  double y = q[2];
  double z = q[3];
  double ww = w * w;
  double xx = x * x;
  double yy = y * y;
  double zz = z * z;
  double xy = x * y;
  double zw = z * w;
  double xz = x * z;
  double yw = y * w;
  double yz = y * z;
  double xw = x * w;

  Q[0 + 0] = ww + xx - yy - zz;
  Q[3 + 0] = 2 * (xy + zw);
  Q[6 + 0] = 2 * (xz - yw);

  Q[0 + 1] = 2 * (xy - zw);
  Q[3 + 1] = ww - xx + yy - zz;
  Q[6 + 1] = 2 * (yz + xw);

  Q[0 + 2] = 2 * (xz + yw);
  Q[3 + 2] = 2 * (yz - xw);
  Q[6 + 2] = ww - xx - yy + zz;
}

void qmat_skewmat(double* S, const double* x) {
  S[0] = 0.0;
  S[1] = +x[2];
  S[2] = -x[1];

  S[3] = -x[2];
  S[4] = 0.0;
  S[5] = +x[0];

  S[6] = +x[1];
  S[7] = -x[0];
  S[8] = 0.0;
}

void qmat_lmat(double* L, const double* q) {
  L[0] = q[0];
  L[1] = q[1];
  L[2] = q[2];
  L[3] = q[3];

  L[4] = -q[1];
  L[5] = q[0];
  L[6] = q[3];
  L[7] = -q[2];

  L[8] = -q[2];
  L[9] = -q[3];
  L[10] = q[0];
  L[11] = q[1];

  L[12] = -q[3];
  L[13] = q[2];
  L[14] = -q[1];
  L[15] = q[0];
}

void qmat_rmat(double* R, const double* q) {
  R[0] = q[0];
  R[1] = q[1];
  R[2] = q[2];
  R[3] = q[3];

  R[4] = -q[1];
  R[5] = q[0];
  R[6] = -q[3];
  R[7] = q[2];

  R[8] = -q[2];
  R[9] = q[3];
  R[10] = q[0];
  R[11] = -q[1];

  R[12] = -q[3];
  R[13] = -q[2];
  R[14] = q[1];
  R[15] = q[0];
}

void qmat_gmat(double* G, double* q) {
  G[0] = -q[1];
  G[1] = q[0];
  G[2] = q[3];
  G[3] = -q[2];

  G[4] = -q[2];
  G[5] = -q[3];
  G[6] = q[0];
  G[7] = q[1];

  G[8] = -q[3];
  G[9] = q[2];
  G[10] = -q[1];
  G[11] = q[0];
}

void qmat_cay(double* q, const double* phi) {
  double m = 1.0 / sqrt(1 + phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
  q[0] = m;
  q[1] = phi[0] * m;
  q[2] = phi[1] * m;
  q[3] = phi[2] * m;
}

void qmat_icay(double* phi, const double* q) {
  phi[0] = q[1] / q[0];
  phi[1] = q[2] / q[0];
  phi[2] = q[3] / q[0];
}

void qmat_dcay(double* D, const double* phi) {
  double m = 1.0 / sqrt(1 + phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
  double m3 = -m * m * m;
  D[0] = m3 * phi[0];
  D[1] = m3 * phi[0] * phi[0];
  D[2] = m3 * phi[1] * phi[0];
  D[3] = m3 * phi[2] * phi[0];

  D[4] = m3 * phi[1];
  D[5] = m3 * phi[0] * phi[1];
  D[6] = m3 * phi[1] * phi[1];
  D[7] = m3 * phi[2] * phi[1];

  D[8] = m3 * phi[2];
  D[9] = m3 * phi[0] * phi[2];
  D[10] = m3 * phi[1] * phi[2];
  D[11] = m3 * phi[2] * phi[2];
  D[1] += m;
  D[6] += m;
  D[11] += m;
}

void qmat_rotate(double* x2, const double* q, const double* x) {
  x2[0] = (q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3]) * x[0];
  x2[0] += 2 * (q[1] * q[2] - q[0] * q[3]) * x[1];
  x2[0] += 2 * (q[0] * q[2] + q[1] * q[3]) * x[2];

  x2[1] = (q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3]) * x[1];
  x2[1] += 2 * (q[1] * q[2] + q[0] * q[3]) * x[0];
  x2[1] += 2 * (q[2] * q[3] - q[0] * q[1]) * x[2];

  x2[2] = (q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]) * x[2];
  x2[2] += 2 * (q[0] * q[1] + q[2] * q[3]) * x[1];
  x2[2] += 2 * (q[1] * q[3] - q[0] * q[2]) * x[0];
}

void star_QuatRotateActiveJacobian(double* D, const double* q, const double* x) {
  D[0] = 2 * q[0] * x[0] + 2 * q[2] * x[2] - 2 * q[3] * x[1];
  D[1] = 2 * q[3] * x[0] + 2 * q[0] * x[1] - 2 * q[1] * x[2];
  D[2] = 2 * q[0] * x[2] + 2 * q[1] * x[1] - 2 * q[2] * x[0];

  D[3] = 2 * q[1] * x[0] + 2 * q[2] * x[1] + 2 * q[3] * x[2];
  D[4] = 2 * q[2] * x[0] - 2 * q[0] * x[2] - 2 * q[1] * x[1];
  D[5] = 2 * q[3] * x[0] + 2 * q[0] * x[1] - 2 * q[1] * x[2];

  D[6] = 2 * q[0] * x[2] + 2 * q[1] * x[1] - 2 * q[2] * x[0];
  D[7] = 2 * q[1] * x[0] + 2 * q[2] * x[1] + 2 * q[3] * x[2];
  D[8] = 2 * q[3] * x[1] - 2 * q[0] * x[0] - 2 * q[2] * x[2];

  D[9] = 2 * q[1] * x[2] - 2 * q[3] * x[0] - 2 * q[0] * x[1];
  D[10] = 2 * q[0] * x[0] + 2 * q[2] * x[2] - 2 * q[3] * x[1];
  D[11] = 2 * q[1] * x[0] + 2 * q[2] * x[1] + 2 * q[3] * x[2];
}

void qmat_err(double* phi, const double* q1, const double* q2) {
  double n = 1 / (q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2] + q1[3] * q2[3]);
  phi[0] = (q1[1] * q2[0] + q1[2] * q2[3] - q1[0] * q2[1] - q1[3] * q2[2]) * n;
  phi[1] = (q1[2] * q2[0] + q1[3] * q2[1] - q1[0] * q2[2] - q1[1] * q2[3]) * n;
  phi[2] = (q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1] - q1[0] * q2[3]) * n;
}
void qmat_adderr(double* q, const double* q0, const double* phi) {
  double dq[4];
  qmat_cay(dq, phi);
  star_QuatCompose(q, q0, dq);
}

void star_QuatToRodriguesParam(double* g, const double* q) {
  double s = q[0];
  if (fabs(s) < STAR_EPS) {
    g[0] = NAN;
    g[1] = NAN;
    g[2] = NAN;
    return;
  }
  g[0] = q[1] / s;
  g[1] = q[2] / s;
  g[2] = q[3] / s;
}

void star_RodriguesParamToQuat(double* q, const double* g) {
  double M = 1.0 / sqrt(1 + g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
  q[0] = M;
  q[1] = g[0] * M;
  q[2] = g[1] * M;
  q[3] = g[2] * M;
}

void star_QuatToMRP(double* p, const double* q) {
  double s = q[0];

  if (fabs(s + 1) < STAR_EPS) {
    p[0] = NAN;
    p[1] = NAN;
    p[2] = NAN;
    return;
  }
  p[0] = q[1] / (1 + s);
  p[1] = q[2] / (1 + s);
  p[2] = q[3] / (1 + s);
}

void star_MRPToQuat(double* q, const double* p) {
  double norm2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  double M = 2.0 / (1 + norm2);
  q[0] = (1 - norm2) / (1 + norm2);
  q[1] = p[0] * M;
  q[2] = p[1] * M;
  q[3] = p[2] * M;
}

void star_QuatToAxisAngle(double* aa, const double* q) {
  star_QuatLogm(aa + 1, q);
  double theta = sqrt(aa[1] * aa[1] + aa[2] * aa[2] + aa[3] * aa[3]);
  if (fabs(theta) < STAR_EPS) {
    aa[0] = 0;
    aa[1] = 1;
    aa[2] = 0;
    aa[3] = 0;
    return;
  }
  aa[0] = theta;
  aa[1] /= theta;
  aa[2] /= theta;
  aa[3] /= theta;
}

void star_QuatToEulerXYZ(double* e, const double* q) {
  double w = q[0];
  double x = q[1];
  double y = q[2];
  double z = q[3];
  double ww = w * w;
  double xx = x * x;
  double yy = y * y;
  double zz = z * z;
  double xy = x * y;
  double zw = z * w;
  double xz = x * z;
  double yw = y * w;
  double yz = y * z;
  double xw = x * w;

  double Q23 = 2 * (yz - xw);
  double Q33 = ww - xx - yy + zz;

  double Q11 = ww + xx - yy - zz;
  double Q12 = 2 * (xy - zw);
  double Q13 = 2 * (xz + yw);

  e[0] = atan2(-Q23, Q33);
  e[1] = atan2(Q13, sqrt(Q11 * Q11 + Q12 * Q12));
  e[2] = atan2(-Q12, Q11);
}

void star_QuatToEulerZYX(double* e, const double* q) {
  double w = q[0];
  double x = q[1];
  double y = q[2];
  double z = q[3];
  double ww = w * w;
  double xx = x * x;
  double yy = y * y;
  double zz = z * z;
  double xy = x * y;
  double zw = z * w;
  double xz = x * z;
  double yw = y * w;
  double yz = y * z;
  double xw = x * w;

  double Q11 = ww + xx - yy - zz;
  double Q21 = 2 * (xy + zw);
  double Q31 = 2 * (xz - yw);
  double Q32 = 2 * (yz + xw);
  double Q33 =  ww - xx - yy + zz;

  e[0] = atan2(Q21, Q11);
  e[1] = atan2(-Q31, sqrt(Q32 * Q32 + Q33 * Q33));
  e[2] = atan2(Q32, Q33);
}
