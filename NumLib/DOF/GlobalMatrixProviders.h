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

#ifdef MSVC
#include "numlib_export.h"
#else
#define NUMLIB_EXPORT
#endif

namespace MathLib
{

template<typename Vector>
struct GlobalVectorProvider
{
	static NUMLIB_EXPORT VectorProvider<Vector>& provider;
};

template<typename Matrix>
struct GlobalMatrixProvider
{
	static NUMLIB_EXPORT MatrixProvider<Matrix>& provider;
};

void cleanupGlobalMatrixProviders();

} // MathLib

#endif // MATHLIB_GLOBAL_MATRIX_PROVIDERS
