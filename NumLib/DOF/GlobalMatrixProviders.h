/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MatrixProviderUser.h"

namespace NumLib
{
struct GlobalVectorProvider
{
    static VectorProvider& provider;
};

struct GlobalMatrixProvider
{
    static MatrixProvider& provider;
};

void cleanupGlobalMatrixProviders();

} // MathLib
