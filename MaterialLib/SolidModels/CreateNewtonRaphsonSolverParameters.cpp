/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 10, 2018, 11:32 AM
 */

#include "CreateNewtonRaphsonSolverParameters.h"

#include "BaseLib/ConfigTree.h"

#include "NumLib/NewtonRaphson.h"

namespace MaterialLib
{
NumLib::NewtonRaphsonSolverParameters createNewtonRaphsonSolverParameters(
    BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const maximum_iterations =
        //! \ogs_file_param{nonlinear_solver__maximum_iterations}
        config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: {:d}.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{nonlinear_solver__error_tolerance}
        config.getConfigParameter<double>("error_tolerance");

    DBUG("\terror_tolerance: {:g}.", error_tolerance);

    return {maximum_iterations, error_tolerance};
}
}  // namespace MaterialLib
