/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateTimeLoop.h"

#include <range/v3/algorithm/any_of.hpp>

#include "BaseLib/ConfigTree.h"
#include "ProcessLib/CreateProcessData.h"
#include "ProcessLib/Output/CreateOutput.h"
#include "ProcessLib/Output/Output.h"
#include "ProcessLib/Output/SubmeshResiduumOutputConfig.h"
#include "TimeLoop.h"

namespace ProcessLib
{
std::unique_ptr<TimeLoop> createTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    const std::vector<std::unique_ptr<Process>>& processes,
    const std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>>&
        nonlinear_solvers,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes,
    bool const compensate_non_equilibrium_initial_residuum)
{
    auto const& coupling_config
        //! \ogs_file_param{prj__time_loop__global_process_coupling}
        = config.getConfigSubtreeOptional("global_process_coupling");

    std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>
        global_coupling_conv_criteria;
    std::map<std::string, int> local_coupling_processes;
    int max_coupling_iterations = 1;
    if (coupling_config)
    {
        max_coupling_iterations
            //! \ogs_file_param{prj__time_loop__global_process_coupling__max_iter}
            = coupling_config->getConfigParameter<int>("max_iter");

        auto const& coupling_convergence_criteria_config =
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criteria}
            coupling_config->getConfigSubtree("convergence_criteria");

        auto coupling_convergence_criterion_config =
            //! \ogs_file_param{prj__time_loop__global_process_coupling__convergence_criteria__convergence_criterion}
            coupling_convergence_criteria_config.getConfigSubtreeList(
                "convergence_criterion");
        std::transform(coupling_convergence_criterion_config.begin(),
                       coupling_convergence_criterion_config.end(),
                       std::back_inserter(global_coupling_conv_criteria),
                       [](BaseLib::ConfigTree const& c)
                       { return NumLib::createConvergenceCriterion(c); });

        auto const& local_coupling_processes_config =
            coupling_config
                //! \ogs_file_param{prj__time_loop__global_process_coupling__local_coupling_processes}
                ->getConfigSubtreeOptional("local_coupling_processes");
        if (local_coupling_processes_config)
        {
            for (
                auto name :
                local_coupling_processes_config
                    //! \ogs_file_param{prj__time_loop__global_process_coupling__local_coupling_processes__process_name}
                    ->getConfigParameterList<std::string>("process_name"))
            {
                BaseLib::insertIfKeyUniqueElseError(
                    local_coupling_processes, name, -1,
                    "The name of locally coupled process is not unique.");
            }

            static std::string const error_info =
                "Please check the number of elements in the tag "
                "'time_loop/global_process_coupling/"
                "local_coupling_processes' in the project file.";
            if (local_coupling_processes.size() >
                global_coupling_conv_criteria.size() - 1)
            {
                OGS_FATAL(
                    "The number of the locally coupled processes is greater "
                    "than or equal to the number of total coupled processes. "
                    "{}",
                    error_info);
            }

            INFO("There are {:d} locally coupled processes.",
                 local_coupling_processes.size());

            if (local_coupling_processes.size() < 2)
            {
                OGS_FATAL(
                    "At least two locally coupled processes are required for "
                    "the local iteration solution! {} ",
                    error_info);
            }
        }
    }

    //! \ogs_file_param{prj__time_loop__output}
    auto output_config_tree = config.getConfigSubtreeOptional("output");
    if (!output_config_tree)
    {
        INFO("No output section found.");
    }
    auto outputs =
        output_config_tree
            ? createOutput(*output_config_tree, output_directory, meshes)
            //! \ogs_file_param{prj__time_loop__outputs}
            : createOutputs(config.getConfigSubtree("outputs"),
                            output_directory, meshes);
    auto const fixed_times_for_output =
        calculateUniqueFixedTimesForAllOutputs(outputs);

    if (auto const submesh_residuum_output_config_tree =
            //! \ogs_file_param{prj__time_loop__submesh_residuum_output}
        config.getConfigSubtreeOptional("submesh_residuum_output");
        submesh_residuum_output_config_tree)
    {
        auto smroc = createSubmeshResiduumOutputConfig(
            *submesh_residuum_output_config_tree, output_directory, meshes);

        for (auto& process : processes)
        {
            auto const& residuum_vector_names =
                process->initializeAssemblyOnSubmeshes(smroc.meshes);

            for (auto const& name : residuum_vector_names)
            {
                smroc.output.doNotProjectFromBulkMeshToSubmeshes(
                    name, MeshLib::MeshItemType::Node);
            }
        }

        outputs.push_back(std::move(smroc.output));
    }
    else
    {
        // Submesh assembly must always be initialized.
        for (auto& process : processes)
        {
            process->initializeAssemblyOnSubmeshes({});
        }
    }

    auto per_process_data = createPerProcessData(
        //! \ogs_file_param{prj__time_loop__processes}
        config.getConfigSubtree("processes"), processes, nonlinear_solvers,
        compensate_non_equilibrium_initial_residuum, fixed_times_for_output,
        local_coupling_processes);

    const bool use_staggered_scheme =
        ranges::any_of(processes.begin(), processes.end(),
                       [](auto const& process)
                       { return !(process->isMonolithicSchemeUsed()); });

    if (!use_staggered_scheme && per_process_data.size() > 1)
    {
        OGS_FATAL(
            "The monolithic scheme is used. However more than one process data "
            "tags (by name \"process\") inside tag \"time_loop\" are defined "
            "for the staggered scheme. If you want to use staggered scheme, "
            "please set the element of tag \"<coupling_scheme>\" to "
            "\"staggered\".");
    }

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
        std::move(outputs), std::move(per_process_data),
        max_coupling_iterations, std::move(global_coupling_conv_criteria),
        std::move(local_coupling_processes), start_time, end_time);
}
}  // namespace ProcessLib
