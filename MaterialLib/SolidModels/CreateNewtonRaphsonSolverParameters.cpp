/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateNewtonRaphsonSolverParameters.cpp
 *  Created on July 10, 2018, 11:32 AM
 */

#include "CreateNewtonRaphsonSolverParameters.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "NumLib/NewtonRaphson.h"

namespace MaterialLib
{
NumLib::NewtonRaphsonSolverParameters createNewtonRaphsonSolverParameters(
    BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");

    auto const maximum_iterations =
        //! \ogs_file_param{material__solid__constitutive_relation__nonlinear_solver__maximum_iterations}
        nonlinear_solver_config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: %d.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{material__solid__constitutive_relation__nonlinear_solver__error_tolerance}
        nonlinear_solver_config.getConfigParameter<double>("error_tolerance");

    DBUG("\terror_tolerance: %g.", error_tolerance);

    return {maximum_iterations, error_tolerance};
}
}
