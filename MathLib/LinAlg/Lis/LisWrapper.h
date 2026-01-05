// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include <lis.h>
#ifdef conj
#undef conj
#endif
