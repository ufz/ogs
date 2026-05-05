// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "NonlinearSolver.h"

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
/*! Creates a new nonlinear solver from the given configuration.
 *
 * \param linear_solver the linear solver that will be used by the nonlinear
 *                      solver
 * \param config configuration settings
 *
 * \return a pair <tt>(nl_slv, tag)</tt> where \c nl_slv is the generated
 *         nonlinear solver instance and the \c tag indicates if it uses
 *         the Picard or Newton-Raphson method
 */
std::pair<std::unique_ptr<NonlinearSolverBase>, NonlinearSolverTag>
createNonlinearSolver(GlobalLinearSolver& linear_solver,
                      BaseLib::ConfigTree const& config);

}  // namespace NumLib
