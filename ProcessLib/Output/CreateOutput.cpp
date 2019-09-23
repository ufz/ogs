/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateOutput.h"

#include <logog/include/logog.hpp>
#include <memory>
#include <tuple>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/Mesh.h"

#include "Output.h"

namespace ProcessLib
{
std::unique_ptr<Output> createOutput(
    const BaseLib::ConfigTree& config,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    DBUG("Parse output configuration:");

    //! \ogs_file_param{prj__time_loop__output__type}
    config.checkConfigParameter("type", "VTK");

    auto const prefix =
        //! \ogs_file_param{prj__time_loop__output__prefix}
        config.getConfigParameter<std::string>("prefix");

    auto const compress_output =
        //! \ogs_file_param{prj__time_loop__output__compress_output}
        config.getConfigParameter("compress_output", true);

    auto const data_mode =
        //! \ogs_file_param{prj__time_loop__output__data_mode}
        config.getConfigParameter<std::string>("data_mode", "Appended");

    // Construction of output times
    std::vector<Output::PairRepeatEachSteps> repeats_each_steps;

    std::vector<double> fixed_output_times;

    //! \ogs_file_param{prj__time_loop__output__timesteps}
    if (auto const timesteps = config.getConfigSubtreeOptional("timesteps"))
    {
        //! \ogs_file_param{prj__time_loop__output__timesteps__pair}
        for (auto pair : timesteps->getConfigSubtreeList("pair"))
        {
            //! \ogs_file_param{prj__time_loop__output__timesteps__pair__repeat}
            auto repeat = pair.getConfigParameter<unsigned>("repeat");
            //! \ogs_file_param{prj__time_loop__output__timesteps__pair__each_steps}
            auto each_steps = pair.getConfigParameter<unsigned>("each_steps");

            assert(repeat != 0 && each_steps != 0);
            repeats_each_steps.emplace_back(repeat, each_steps);
        }

        if (repeats_each_steps.empty())
        {
            OGS_FATAL(
                "You have not given any pair (<repeat/>, <each_steps/>) that "
                "defines"
                " at which timesteps output shall be written. Aborting.");
        }
    }
    else
    {
        repeats_each_steps.emplace_back(1, 1);
    }

    //! \ogs_file_param{prj__time_loop__output__variables}
    auto const out_vars = config.getConfigSubtree("variables");

    std::set<std::string> output_variables;
    for (auto out_var :
         //! \ogs_file_param{prj__time_loop__output__variables__variable}
         out_vars.getConfigParameterList<std::string>("variable"))
    {
        if (output_variables.find(out_var) != output_variables.cend())
        {
            OGS_FATAL("output variable `%s' specified more than once.",
                      out_var.c_str());
        }

        DBUG("adding output variable `%s'", out_var.c_str());
        output_variables.insert(out_var);
    }

    //! \ogs_file_param{prj__time_loop__output__output_extrapolation_residuals}
    bool const output_residuals = config.getConfigParameter<bool>(
        "output_extrapolation_residuals", false);

    ProcessOutput process_output{output_variables, output_residuals};

    std::vector<std::string> mesh_names_for_output;
    //! \ogs_file_param{prj__time_loop__output__meshes}
    if (auto const meshes_config = config.getConfigSubtreeOptional("meshes"))
    {
        //! \ogs_file_param{prj__time_loop__output__meshes__mesh}
        for (auto mesh_config : meshes_config->getConfigParameterList("mesh"))
        {
            mesh_names_for_output.push_back(
                mesh_config.getValue<std::string>());
            INFO("Configure mesh '%s' for output.",
                 mesh_names_for_output.back().c_str());
        }
    }

    auto fixed_output_times_ptr =
        //! \ogs_file_param{prj__time_loop__output__fixed_output_times}
        config.getConfigParameterOptional<std::vector<double>>(
            "fixed_output_times");
    if (fixed_output_times_ptr)
    {
        fixed_output_times = std::move(*fixed_output_times_ptr);
        // Remove possible duplicated elements and sort in descending order.
        BaseLib::makeVectorUnique(fixed_output_times, std::greater<>());
    }

    bool const output_iteration_results =
        //! \ogs_file_param{prj__time_loop__output__output_iteration_results}
        config.getConfigParameter<bool>("output_iteration_results", false);

    return std::make_unique<Output>(
        output_directory, prefix, compress_output, data_mode,
        output_iteration_results, std::move(repeats_each_steps),
        std::move(fixed_output_times), std::move(process_output),
        std::move(mesh_names_for_output), meshes);
}

}  // namespace ProcessLib
