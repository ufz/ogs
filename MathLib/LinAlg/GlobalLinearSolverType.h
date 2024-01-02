/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#if defined(USE_LIS)

#include "MathLib/LinAlg/EigenLis/EigenLisLinearSolver.h"

using GlobalLinearSolver = MathLib::EigenLisLinearSolver;

#elif defined(USE_PETSC)

#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"

using GlobalLinearSolver = MathLib::PETScLinearSolver;

#else

#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

using GlobalLinearSolver = MathLib::EigenLinearSolver;

#endif
