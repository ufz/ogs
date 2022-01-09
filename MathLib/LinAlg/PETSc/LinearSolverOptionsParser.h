/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/LinearSolverOptions.h"
#include "MathLib/LinAlg/LinearSolverOptionsParser.h"
#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"

namespace MathLib
{
template <>
struct LinearSolverOptionsParser<PETScLinearSolver> final
{
    /// The method parses the linear solver options for PETSc.
    /// @param solver_prefix the prefix for the linear solver to distinguish
    /// different linear solvers for instance in the staggered schema
    /// @param config the part of the property tree (usually created from the
    /// linear solver section in the project file)
    /// @return the first item of the returned tuple is the solver prefix as
    /// string, the second item are all the options passed as string to PETSc
    /// options database
    std::tuple<std::string, std::string> parseNameAndOptions(
        std::string solver_prefix,
        BaseLib::ConfigTree const* const config) const
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
                    subtree->getConfigParameterOptional<std::string>(
                        "parameters"))
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
};

}  // namespace MathLib
