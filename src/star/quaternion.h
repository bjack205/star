//
// Created by Brian Jackson on 3/12/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

void qmat_inv(double* qinv, const double* q);
void qmat_quat2rotmat(double* Q, const double* q);
void qmat_skewmat(double* S, const double* x);
void qmat_lmat(double* L, const double* q);
void qmat_rmat(double* R, const double* q);
void qmat_gmat(double* G, double* q);
void qmat_cay(double* q, const double* phi);
void qmat_icay(double* phi, const double* q);
void qmat_dcay(double* D, const double* phi);
void qmat_rotate(double* x2, const double* q, const double* x);
void qmat_drotate(double* D, const double* q, const double* x);
void qmat_compose(double* q3, const double* q1, const double* q2);
void qmat_err(double* phi, const double* q1, const double* q2);
void qmat_adderr(double* q, const double* q0, const double* phi);
