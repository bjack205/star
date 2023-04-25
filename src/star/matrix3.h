//
// Created by Brian Jackson on 3/9/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "typedefs.h"

/*---------------------------------*/
/* Setters                         */
/*---------------------------------*/

void star_SetZero33(sfloat mat[9]);
void star_SetConst33(sfloat mat[9], sfloat value);
void star_SetIdentity33(sfloat mat[9], sfloat value);
void star_SetDiagonal33(sfloat mat[9], const sfloat diag[3]);
void star_Copy33(sfloat dst[9], const sfloat src[9]);
void star_Transpose33(sfloat dst[9], const sfloat src[9]);
void star_TransposeInPlace33(sfloat mat[9]);

/*---------------------------------*/
/* Multiplication                  */
/*---------------------------------*/

void star_MatMul33(sfloat C[9], const sfloat A[9], const sfloat B[9]);
void star_VecMul33(sfloat C[3], const sfloat A[9], const sfloat x[3]);

void star_TransposedMatMul33(sfloat C[9], const sfloat At[9], const sfloat B[9]);
void star_MatMulTransposed33(sfloat C[9], const sfloat A[9], const sfloat Bt[9]);
void star_TransposedVecMul33(sfloat y[3], const sfloat At[9], const sfloat x[3]);

/*---------------------------------*/
/* Triangular Matrices             */
/*---------------------------------*/

void star_UpperMatMul33(sfloat C[9], const sfloat U[9], const sfloat A[9]);
void star_UpperVecMul33(sfloat C[3], const sfloat U[9], const sfloat x[3]);
void star_LowerMatMul33(sfloat C[9], const sfloat L[9], const sfloat A[9]);
void star_LowerVecMul33(sfloat C[3], const sfloat L[9], const sfloat x[3]);

void star_UpperTriSolve33(sfloat x[3], const sfloat U[9], const sfloat b[3]);
void star_LowerTriSolve33(sfloat x[3], const sfloat L[9], const sfloat b[3]);

/*---------------------------------*/
/* Linear Algebra                  */
/*---------------------------------*/

sfloat star_Det33(const sfloat mat[9]);

// Decompositions
void star_Chol33(sfloat U[9], const sfloat mat[9]);
void star_QR33(sfloat Q[9], sfloat R[9], const sfloat A[9]);
void star_LU33(sfloat Q[9], sfloat R[9], const sfloat A[9]);
void star_Eigen33(sfloat eigenvalues[3], sfloat eigenvectors[9], const sfloat mat[9]);
void star_SVD33(sfloat U[9], sfloat S[3], sfloat V[9], const sfloat A[9]);

// Solves
void star_CholSolve33(sfloat x[3], const sfloat A[9], const sfloat b[3]);

// Inverses
void star_Inverse33(sfloat mat[9]);
void star_InversePSD33(sfloat mat[9]);

