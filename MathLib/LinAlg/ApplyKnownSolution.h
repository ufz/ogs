/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef MATHLIB_APPLYKNOWNSOLUTION_H_
#define MATHLIB_APPLYKNOWNSOLUTION_H_

#include <vector>

#ifdef OGS_USE_EIGEN
#include "MathLib/LinAlg/Eigen/EigenTools.h"
#endif // OGS_USE_EIGEN

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScTools.h"
#endif // USE_PETSC

#endif  // MATHLIB_APPLYKNOWNSOLUTION_H_
