/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#ifdef OGS_USE_EIGEN
#include "MathLib/LinAlg/Eigen/EigenTools.h"
#endif // OGS_USE_EIGEN

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScTools.h"
#endif // USE_PETSC
