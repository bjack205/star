//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "typedefs.h"

/*---------------------------------*/
/* Setters                         */
/*---------------------------------*/

void star_SetZero44(sfloat mat[16]);
void star_SetConst44(sfloat mat[16], sfloat value);
void star_SetIdentity44(sfloat mat[16], sfloat val);
void star_SetDiagonal44(sfloat mat[16], const sfloat diag[4]);
void star_Copy44(sfloat dst[16], const sfloat src[16]);
void star_Transpose44(sfloat dst[16], const sfloat src[16]);
void star_TransposeInPlace44(sfloat mat[16]);

/*---------------------------------*/
/* Multiplication                  */
/*---------------------------------*/

void star_MatMul44(sfloat C[16], const sfloat A[16], const sfloat B[16]);
void star_VecMul44(sfloat y[4], const sfloat A[16], const sfloat x[4]);

void star_TransposedMatMul44(sfloat C[16], const sfloat At[16], const sfloat B[16]);
void star_MatMulTransposed44(sfloat C[16], const sfloat A[16], const sfloat Bt[16]);
void star_TransposedVecMul44(sfloat y[4], const sfloat At[16], const sfloat x[4]);

/*---------------------------------*/
/* Element-wise Operations         */
/*---------------------------------*/

void star_Add44(sfloat C[16], const sfloat A[16], const sfloat B[16]);
void star_Sub44(sfloat C[16], const sfloat A[16], const sfloat B[16]);
void star_Mul44(sfloat C[16], const sfloat A[16], const sfloat B[16]);
void star_Div44(sfloat C[16], const sfloat A[16], const sfloat B[16]);

void star_AddConst44(sfloat C[16], const sfloat A[16], sfloat b);
void star_SubConst44(sfloat C[16], const sfloat A[16], sfloat b);
void star_MulConst44(sfloat C[16], const sfloat A[16], sfloat b);
void star_DivConst44(sfloat C[16], const sfloat A[16], sfloat b);
