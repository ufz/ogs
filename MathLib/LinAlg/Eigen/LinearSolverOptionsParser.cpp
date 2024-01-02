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
#include "BaseLib/Error.h"
#include "MathLib/LinAlg/LinearSolverOptions.h"

namespace MathLib
{
std::tuple<std::string, EigenOption>
LinearSolverOptionsParser<EigenLinearSolver>::parseNameAndOptions(
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
        options.solver_type = MathLib::EigenOption::getSolverType(*solver_type);
    }
    if (auto precon_type =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__precon_type}
        config->getConfigParameterOptional<std::string>("precon_type"))
    {
        options.precon_type = MathLib::EigenOption::getPreconType(*precon_type);
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
    if (auto l =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__l}
        config->getConfigParameterOptional<int>("l"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        options.l = *l;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules.");
#endif
    }
    if (auto s =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__s}
        config->getConfigParameterOptional<int>("s"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        options.s = *s;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules.");
#endif
    }
    if (auto smoothing =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__smoothing}
        config->getConfigParameterOptional<int>("smoothing"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        options.smoothing = *smoothing;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules.");
#endif
    }
    if (auto angle =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__angle}
        config->getConfigParameterOptional<int>("angle"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        options.angle = *angle;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules.");
#endif
    }
    if (auto residualupdate =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__residual_update}
        config->getConfigParameterOptional<int>("residual_update"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        options.residualupdate = *residualupdate;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules.");
#endif
    }
    return {prefix, options};
}
}  // namespace MathLib
