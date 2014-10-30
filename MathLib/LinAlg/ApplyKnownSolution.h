/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef MATHLIB_APPLYKNOWNSOLUTION_H_
#define MATHLIB_APPLYKNOWNSOLUTION_H_

#include <vector>

#include "MathLib/LinAlg/Dense/DenseTools.h"

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisTools.h"
#endif // USE_LIS

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScTools.h"
#endif // USE_PETSC

#endif  // MATHLIB_APPLYKNOWNSOLUTION_H_
