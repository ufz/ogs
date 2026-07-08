// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateNonlinearSolver.h"

#include <spdlog/fmt/ranges.h>

#include <boost/algorithm/string.hpp>
#include <memory>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "DampingReductionStrategy.h"
#include "FixedDampingStrategy.h"
#include "PETScNonlinearSolver.h"

namespace NumLib
{
std::pair<std::unique_ptr<NonlinearSolverBase>, NonlinearSolverTag>
createNonlinearSolver(GlobalLinearSolver& linear_solver,
                      BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__type}
    auto const type = config.getConfigParameter<std::string>("type");

    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__max_iter}
    auto const max_iter = config.getConfigParameter<int>("max_iter");

    if (type == "Picard")
    {
        //! \ogs_file_param_special{prj__nonlinear_solvers__nonlinear_solver__Picard}
        auto const damping =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__Picard__damping}
            config.getConfigParameter<double>("damping", 1.0);
        if (damping <= 0.0 || damping > 1.0)
        {
            OGS_FATAL(
                "The damping factor for the Picard method must be in (0, 1], "
                "got {:g}.",
                damping);
        }
        auto const tag = NonlinearSolverTag::Picard;
        using ConcreteNLS = NonlinearSolver<tag>;
        return std::make_pair(
            std::make_unique<ConcreteNLS>(linear_solver, max_iter, damping),
            tag);
    }
    if (type == "Newton")
    {
        auto const recompute_jacobian =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__Newton__recompute_jacobian}
            config.getConfigParameter<int>("recompute_jacobian", 1);

        //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__Newton__damping}
        auto const damping = config.getConfigParameter<double>("damping", 1.0);
        if (damping <= 0)
        {
            OGS_FATAL(
                "The damping factor for the Newton method must be positive, "
                "got "
                "{:g}.",
                damping);
        }
        auto const damping_reduction =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__Newton__damping_reduction}
            config.getConfigParameterOptional<double>("damping_reduction");
        std::unique_ptr<NewtonStepStrategy> standard_newton;
        if (damping_reduction)
        {
            standard_newton = std::make_unique<DampingReductionStrategy>(
                damping, *damping_reduction);
        }
        else
        {
            standard_newton = std::make_unique<FixedDampingStrategy>(damping);
        }
        auto const tag = NonlinearSolverTag::Newton;
        using ConcreteNLS = NonlinearSolver<tag>;
        auto solver = std::make_unique<ConcreteNLS>(linear_solver, max_iter,
                                                    std::move(standard_newton),
                                                    recompute_jacobian);
        auto const tikhonov_config =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__tikhonov}
            config.getConfigSubtreeOptional("tikhonov");
        if (tikhonov_config)
        {
            auto const tikhonov_lambda =
                //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__tikhonov__lambda}
                tikhonov_config->getConfigParameter<double>("lambda");
            auto const tikhonov_starting_iteration =
                //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__tikhonov__starting_iteration}
                tikhonov_config->getConfigParameter<int>("starting_iteration",
                                                         0);
            solver->setTikhonovLambda(tikhonov_lambda,
                                      tikhonov_starting_iteration);
            DBUG("Newton solver: Tikhonov regularization lambda = {:g}",
                 tikhonov_lambda);
        }
        return std::make_pair(std::move(solver), tag);
    }
#ifdef USE_PETSC
    if (boost::iequals(type, "PETScSNES"))
    {
        auto prefix =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__PETScSNES__prefix}
            config.getConfigParameter<std::string>("prefix", "");
        auto const tag = NonlinearSolverTag::Newton;
        using ConcreteNLS = PETScNonlinearSolver;
        return std::make_pair(std::make_unique<ConcreteNLS>(
                                  linear_solver, max_iter, std::move(prefix)),
                              tag);
    }
    static constexpr std::array valid_types = {"PETScSNES", "Newton", "Picard"};
#else
    static constexpr std::array valid_types = {"Newton", "Picard"};
#endif

    OGS_FATAL("Invalid non-linear solver type '{}'. Supported values: {}.",
              type, fmt::join(valid_types, ", "));
}

}  // namespace NumLib
