/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/Algorithm.h"

#include "NumLib/ODESolver/TimeDiscretizationBuilder.h"
#include "NumLib/TimeStepping/CreateTimeStepper.h"

#include "CreateProcessData.h"

namespace ProcessLib
{
static std::unique_ptr<ProcessData> makeProcessData(
    std::unique_ptr<NumLib::TimeStepAlgorithm>&& timestepper,
    NumLib::NonlinearSolverBase& nonlinear_solver,
    int const process_id,
    Process& process,
    std::unique_ptr<NumLib::TimeDiscretization>&& time_disc,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& conv_crit,
    bool const compensate_initial_out_of_balance_forces)
{
    using Tag = NumLib::NonlinearSolverTag;

    if (auto* nonlinear_solver_picard =
            dynamic_cast<NumLib::NonlinearSolver<Tag::Picard>*>(
                &nonlinear_solver))
    {
        nonlinear_solver_picard->compensateInitialOutOfBalanceForces(
            compensate_initial_out_of_balance_forces);
        return std::make_unique<ProcessData>(
            std::move(timestepper), *nonlinear_solver_picard,
            std::move(conv_crit), std::move(time_disc), process_id, process);
    }
    if (auto* nonlinear_solver_newton =
            dynamic_cast<NumLib::NonlinearSolver<Tag::Newton>*>(
                &nonlinear_solver))
    {
        nonlinear_solver_newton->compensateInitialOutOfBalanceForces(
            compensate_initial_out_of_balance_forces);
        return std::make_unique<ProcessData>(
            std::move(timestepper), *nonlinear_solver_newton,
            std::move(conv_crit), std::move(time_disc), process_id, process);
    }

    OGS_FATAL("Encountered unknown nonlinear solver type. Aborting");
}

std::vector<std::unique_ptr<ProcessData>> createPerProcessData(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<Process>> const& processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers)
{
    std::vector<std::unique_ptr<ProcessData>> per_process_data;
    int process_id = 0;

    //! \ogs_file_param{prj__time_loop__processes__process}
    for (auto pcs_config : config.getConfigSubtreeList("process"))
    {
        //! \ogs_file_attr{prj__time_loop__processes__process__ref}
        auto const pcs_name = pcs_config.getConfigAttribute<std::string>("ref");
        auto& pcs = *BaseLib::getIfOrError(
            processes,
            [&pcs_name](std::unique_ptr<Process> const& p) {
                return p->name == pcs_name;
            },
            "A process with the given name has not been defined.");

        auto const nl_slv_name =
            //! \ogs_file_param{prj__time_loop__processes__process__nonlinear_solver}
            pcs_config.getConfigParameter<std::string>("nonlinear_solver");
        auto& nl_slv = *BaseLib::getOrError(
            nonlinear_solvers, nl_slv_name,
            "A nonlinear solver with the given name has not been defined.");

        auto time_disc = NumLib::createTimeDiscretization(
            //! \ogs_file_param{prj__time_loop__processes__process__time_discretization}
            pcs_config.getConfigSubtree("time_discretization"));

        auto timestepper = NumLib::createTimeStepper(
            //! \ogs_file_param{prj__time_loop__processes__process__time_stepping}
            pcs_config.getConfigSubtree("time_stepping"));

        auto conv_crit = NumLib::createConvergenceCriterion(
            //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion}
            pcs_config.getConfigSubtree("convergence_criterion"));

        auto const compensate_initial_out_of_balance_forces =
            //! \ogs_file_param{prj__time_loop__processes__process__compensate_initial_out_of_balance_forces}
            pcs_config.getConfigParameter<bool>(
                "compensate_initial_out_of_balance_forces", false);

        //! \ogs_file_param{prj__time_loop__processes__process__output}
        auto output = pcs_config.getConfigSubtreeOptional("output");
        if (output)
        {
            OGS_FATAL(
                "In order to make the specification of output in the project "
                "file consistent, the variables output tags were moved from "
                "xpath "
                "'//OpenGeoSysProject/time_loop/processes/process/output' "
                "to the global output section, i.e., to the xpath "
                "'//OpenGeoSysProject/time_loop/output'. This has to be done "
                "in the current project file!");
        }

        per_process_data.emplace_back(
            makeProcessData(std::move(timestepper), nl_slv, process_id, pcs,
                            std::move(time_disc), std::move(conv_crit),
                            compensate_initial_out_of_balance_forces));
        ++process_id;
    }

    if (per_process_data.size() != processes.size())
    {
        if (processes.size() > 1)
        {
            OGS_FATAL(
                "Some processes have not been configured to be solved by this "
                " time loop.");
        }
        else
        {
            INFO(
                "The equations of the coupled processes will be solved by the "
                "staggered scheme.");
        }
    }

    return per_process_data;
}
}  // namespace ProcessLib
