//
// Created by Brian Jackson on 4/22/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#include "typedefs.h"

/*---------------------------------*/
/* Setters                         */
/*---------------------------------*/

void star_SetZero43(sfloat mat[12]);
void star_SetConst43(sfloat mat[12], sfloat value);
void star_Copy43(sfloat dst[12], const sfloat src[12]);
//void star_Transpose43(sfloat dst[12], const sfloat src[12]);

/*---------------------------------*/
/* Multiplication                  */
/*---------------------------------*/

/*
 * @brief Multiply a 4x3 matrix by a 3x3 matrix
 */
void star_MatMul433(sfloat C43[12], const sfloat A43[12], const sfloat B33[9]);
void star_MatMulTransposed433(sfloat C43[12], const sfloat A43[12], const sfloat B33t[9]);

/*
 * @brief Multiply a 4x4 matrix by a 4x3 matrix
 */
void star_MatMul443(sfloat C43[12], const sfloat A44[16], const sfloat B43[12]);
void star_TransposedMatMul443(sfloat C43[12], const sfloat A44t[16], const sfloat B43[12]);

/*
 * @brief Multiply a 3x4 matrix by a 4x4 matrix
 */
void star_MatMul344(sfloat C34[12], const sfloat A34[12], const sfloat B44[16]);
void star_MatMulTransposed344(sfloat C34[12], const sfloat A34[12], const sfloat B44t[16]);

/*
 * @brief Multiply a 3x3 matrix by a 3x4 matrix
 */
void star_MatMul334(sfloat C34[12], const sfloat A33[9], const sfloat B34[12]);
void star_TransposedMatMul334(sfloat C34[12], const sfloat A33t[9], const sfloat B34[12]);

void star_VecMul43(sfloat y[4], const sfloat A[12], const sfloat x[3]);
void star_TransposedVecMul43(sfloat y[3], const sfloat At[12], const sfloat x[4]);





