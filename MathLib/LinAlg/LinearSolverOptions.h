/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include "BaseLib/ConfigTree-fwd.h"

namespace MathLib
{

/*! Ignore linear solver settings not needed for the selected one.
 *
 * The project files support specifying linear solver options for all
 * known solver libraries (currently PETSC, LIS, Eigen) even though for
 * a specific build only one of those settings is used.
 * That clearly conflicts with the requirement of the config tree that
 * each setting present in the project file must be read exactly once.
 *
 * The purpose of this function is to explicitly ignore all the settings
 * that are not relevant for the currently used linear solver
 *
 * \param config The config tree snippet for the linear solver.
 * \param solver_name The tag under which the relevant configuration is found.
 *                    All other configurations will be ignored.
 *
 * This function is currently used in the option parsing code of our
 * \c EigenLinearSolver, \c LisOption and \c PETScLinearSolver
 */
void ignoreOtherLinearSolvers(BaseLib::ConfigTree const& config,
                              std::string const& solver_name);

}  // namespace MathLib
