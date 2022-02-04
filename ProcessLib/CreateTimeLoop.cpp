/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateTimeLoop.h"

#include "BaseLib/ConfigTree.h"
#include "ProcessLib/CreateProcessData.h"
#include "ProcessLib/Output/CreateOutput.h"
#include "ProcessLib/Output/Output.h"
#include "TimeLoop.h"

namespace ProcessLib
{
std::unique_ptr<TimeLoop> createTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    const std::vector<std::unique_ptr<Process>>& processes,
    const std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>>&
        nonlinear_solvers,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    bool const compensate_non_equilibrium_initial_residuum)
{
    auto const& coupling_config
        //! \ogs_file_param{prj__time_loop__global_process_coupling}
        = config.getConfigSubtreeOptional("global_process_coupling");

    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>
        global_coupling_conv_criteria;
    int max_coupling_iterations = 1;
    if (coupling_config)
    {
        max_coupling_iterations
            //! \ogs_file_param{prj__time_loop__global_process_coupling__max_iter}
            = coupling_config->getConfigParameter<int>("max_iter");

        auto const& coupling_convergence_criteria_config =
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criteria}
            coupling_config->getConfigSubtree("convergence_criteria");

        for (
            auto coupling_convergence_criterion_config :
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criteria__convergence_criterion}
            coupling_convergence_criteria_config.getConfigSubtreeList(
                "convergence_criterion"))
        {
            global_coupling_conv_criteria.push_back(
                NumLib::createConvergenceCriterion(
                    coupling_convergence_criterion_config));
        }
    }

    auto output =
        //! \ogs_file_param{prj__time_loop__output}
        createOutput(config.getConfigSubtree("output"), output_directory,
                     meshes);

    auto per_process_data = createPerProcessData(
        //! \ogs_file_param{prj__time_loop__processes}
        config.getConfigSubtree("processes"), processes, nonlinear_solvers,
        compensate_non_equilibrium_initial_residuum);

    if (coupling_config)
    {
        if (global_coupling_conv_criteria.size() != per_process_data.size())
        {
            OGS_FATAL(
                "The number of convergence criteria of the global staggered "
                "coupling loop is not identical to the number of the "
                "processes! Please check the element by tag "
                "global_process_coupling in the project file.");
        }
    }

    const auto minmax_iter =
        std::minmax_element(per_process_data.begin(),
                            per_process_data.end(),
                            [](std::unique_ptr<ProcessData> const& a,
                               std::unique_ptr<ProcessData> const& b) {
                                return (a->timestep_algorithm->end() <
                                        b->timestep_algorithm->end());
                            });
    const double start_time =
        per_process_data[minmax_iter.first - per_process_data.begin()]
            ->timestep_algorithm->begin();
    const double end_time =
        per_process_data[minmax_iter.second - per_process_data.begin()]
            ->timestep_algorithm->end();

    return std::make_unique<TimeLoop>(
        std::move(output), std::move(per_process_data), max_coupling_iterations,
        std::move(global_coupling_conv_criteria), start_time, end_time);
}
}  // namespace ProcessLib
