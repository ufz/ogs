/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "numlib_export.h"

#include "MatrixProviderUser.h"

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

} // MathLib
