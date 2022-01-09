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

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "EigenOption.h"
#include "MathLib/LinAlg/LinearSolverOptions.h"
#include "MathLib/LinAlg/LinearSolverOptionsParser.h"

namespace MathLib
{
template <>
struct LinearSolverOptionsParser<EigenLinearSolver> final
{
    /// The method parses the linear solver options for Eigen.
    /// @param prefix the prefix for the linear solver to distinguish different
    /// linear solvers for instance in the staggered schema
    /// @param solver_config the part of the property tree (usually created from
    /// the linear solver section in the project file)
    /// @return the first item of the returned tuple is the solver prefix as
    /// string, the second item are all the options as an EigenOption object
    std::tuple<std::string, EigenOption> parseNameAndOptions(
        std::string const& prefix,
        BaseLib::ConfigTree const* const solver_config) const
    {
        if (!solver_config)
        {
            return {prefix, EigenOption{}};
        }

        ignoreOtherLinearSolvers(*solver_config, "eigen");
        //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen}
        auto const config = solver_config->getConfigSubtreeOptional("eigen");
        if (!config)
        {
            OGS_FATAL(
                "OGS was compiled with Eigen but the config section in the "
                "project file seems to be invalid");
        }

        EigenOption options;

        if (auto solver_type =
                //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__solver_type}
            config->getConfigParameterOptional<std::string>("solver_type"))
        {
            options.solver_type =
                MathLib::EigenOption::getSolverType(*solver_type);
        }
        if (auto precon_type =
                //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__precon_type}
            config->getConfigParameterOptional<std::string>("precon_type"))
        {
            options.precon_type =
                MathLib::EigenOption::getPreconType(*precon_type);
        }
        if (auto error_tolerance =
                //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__error_tolerance}
            config->getConfigParameterOptional<double>("error_tolerance"))
        {
            options.error_tolerance = *error_tolerance;
        }
        if (auto max_iteration_step =
                //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__max_iteration_step}
            config->getConfigParameterOptional<int>("max_iteration_step"))
        {
            options.max_iterations = *max_iteration_step;
        }
        if (auto scaling =
                //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__scaling}
            config->getConfigParameterOptional<bool>("scaling"))
        {
#ifdef USE_EIGEN_UNSUPPORTED
            options.scaling = *scaling;
#else
            OGS_FATAL(
                "The code is not compiled with the Eigen unsupported modules. "
                "scaling is not available.");
#endif
        }
        if (auto restart =
                //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__restart}
            config->getConfigParameterOptional<int>("restart"))
        {
#ifdef USE_EIGEN_UNSUPPORTED
            options.restart = *restart;
#else
            OGS_FATAL(
                "The code is not compiled with the Eigen unsupported modules. "
                "GMRES/GMRES option restart is not available.");
#endif
        }
        return {prefix, options};
    }
};

}  // namespace MathLib
