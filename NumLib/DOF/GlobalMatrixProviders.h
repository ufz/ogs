/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_GLOBAL_MATRIX_PROVIDERS
#define MATHLIB_GLOBAL_MATRIX_PROVIDERS

#include "MatrixProviderUser.h"

namespace MathLib
{

struct GlobalVectorProvider
{
    static VectorProvider<GlobalVector>& provider;
};

struct GlobalMatrixProvider
{
    static MatrixProvider<GlobalMatrix>& provider;
};

void cleanupGlobalMatrixProviders();

} // MathLib

#endif // MATHLIB_GLOBAL_MATRIX_PROVIDERS
