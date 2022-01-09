/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateNewtonRaphsonSolverParameters.h"

#include "BaseLib/ConfigTree.h"
#include "NewtonRaphson.h"

namespace NumLib
{
NewtonRaphsonSolverParameters createNewtonRaphsonSolverParameters(
    BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const maximum_iterations =
        //! \ogs_file_param{nonlinear_solver__maximum_iterations}
        config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: {:d}.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{nonlinear_solver__error_tolerance}
        config.getConfigParameterOptional<double>("error_tolerance");
    if (error_tolerance)
    {
        WARN(
            "The 'error_tolerance' tag for the Newton-Raphson solver is "
            "deprecated.\n"
            "Use new tags 'residuum_tolerance' and 'increment_tolerance'.\n"
            "For now we use residuum_tolerance {} and increment_tolerance 0.",
            *error_tolerance);
        return {maximum_iterations, *error_tolerance, 0};
    }

    auto const residuum_tolerance =
        //! \ogs_file_param{nonlinear_solver__residuum_tolerance}
        config.getConfigParameter<double>("residuum_tolerance");

    DBUG("\tresiduum_tolerance: {:g}.", residuum_tolerance);

    auto const increment_tolerance =
        //! \ogs_file_param{nonlinear_solver__increment_tolerance}
        config.getConfigParameter<double>("increment_tolerance");

    DBUG("\tincrement_tolerance: {:g}.", increment_tolerance);

    return {maximum_iterations, residuum_tolerance, increment_tolerance};
}
}  // namespace NumLib
