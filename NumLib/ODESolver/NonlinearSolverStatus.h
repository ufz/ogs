/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace NumLib
{
/// Status of the non-linear solver.
struct NonlinearSolverStatus
{
    bool error_norms_met = false;
    int number_iterations = -1;
};
}  // namespace NumLib
