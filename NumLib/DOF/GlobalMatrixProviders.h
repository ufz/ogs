// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MatrixProvider.h"
#include "VectorProvider.h"
#include "numlib_export.h"

namespace NumLib
{
struct GlobalVectorProvider
{
    static NUMLIB_EXPORT VectorProvider& provider;
};

struct GlobalMatrixProvider
{
    static NUMLIB_EXPORT MatrixProvider& provider;
};

void cleanupGlobalMatrixProviders();

}  // namespace NumLib
