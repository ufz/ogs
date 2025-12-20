// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
