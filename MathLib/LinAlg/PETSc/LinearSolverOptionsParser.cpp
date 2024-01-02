/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearSolverOptionsParser.h"

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/LinearSolverOptions.h"

namespace MathLib
{
std::tuple<std::string, std::string>
LinearSolverOptionsParser<PETScLinearSolver>::parseNameAndOptions(
    std::string solver_prefix, BaseLib::ConfigTree const* const config) const
{
    // Insert options into petsc database. Default options are given in the
    // string below.
    std::string petsc_options =
        "-ksp_type cg -pc_type bjacobi -ksp_rtol 1e-16 -ksp_max_it 10000";

    if (config)
    {
        ignoreOtherLinearSolvers(*config, "petsc");

        //! \ogs_file_param{prj__linear_solvers__linear_solver__petsc}
        if (auto const subtree = config->getConfigSubtreeOptional("petsc"))
        {
            if (auto const parameters =
                    //! \ogs_file_param{prj__linear_solvers__linear_solver__petsc__parameters}
                subtree->getConfigParameterOptional<std::string>("parameters"))
            {
                petsc_options = *parameters;
            }

            if (auto const pre =
                    //! \ogs_file_param{prj__linear_solvers__linear_solver__petsc__prefix}
                subtree->getConfigParameterOptional<std::string>("prefix"))
            {
                if (!pre->empty())
                    solver_prefix = *pre + "_";
            }
        }
    }
    return {solver_prefix, petsc_options};
}
}  // namespace MathLib
