/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateTimeLoop.h"

#include <algorithm>
#include <range/v3/algorithm/any_of.hpp>
#include <string>

#include "BaseLib/ConfigTree.h"
#include "NumLib/StaggeredCoupling/CreateStaggeredCoupling.h"
#include "NumLib/StaggeredCoupling/StaggeredCoupling.h"
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
        compensate_non_equilibrium_initial_residuum, fixed_times_for_output);

    const bool use_staggered_scheme =
        ranges::any_of(processes.begin(), processes.end(),
                       [](auto const& process)
                       { return !(process->isMonolithicSchemeUsed()); });

    std::unique_ptr<NumLib::StaggeredCoupling> staggered_coupling = nullptr;
    if (use_staggered_scheme)
    {
        staggered_coupling = NumLib::createStaggeredCoupling<ProcessData>(
            config, per_process_data);
    }
    else
    {
        if (per_process_data.size() > 1)
        {
            OGS_FATAL(
                "The monolithic scheme is used. However more than one "
                "process data tags (by name \"process\") inside tag "
                "\"time_loop\" are defined for the staggered scheme. If you "
                "want to use staggered scheme, please set the element of tag "
                "\"<coupling_scheme>\" to \"staggered\".");
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
        std::move(staggered_coupling), start_time, end_time);
}
}  // namespace ProcessLib
