//
// Created by Brian Jackson on 3/9/23.
// Copyright (c) 2023. All rights reserved.
//

#pragma once

#ifdef STAR_FLOAT
typedef STAR_FLOAT sfloat;
#else
typedef double sfloat;
#endif

#define STAR_EPS 1e-8

// check if c++
#ifdef __cplusplus
#include <utility>
using IndexPair = std::pair<int, int>;
#endif
