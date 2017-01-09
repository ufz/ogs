/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_GLOBAL_MATRIX_PROVIDERS
#define NUMLIB_GLOBAL_MATRIX_PROVIDERS

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

#endif // NUMLIB_GLOBAL_MATRIX_PROVIDERS
