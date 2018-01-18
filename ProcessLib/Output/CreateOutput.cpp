/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateOutput.h"

#include <logog/include/logog.hpp>
#include <memory>
#include <tuple>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"

#include "Output.h"

namespace ProcessLib
{
std::unique_ptr<Output> createOutput(const BaseLib::ConfigTree& config,
                                     std::string const& output_directory)
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
        config.getConfigParameter<std::string>("data_mode", "Binary");

    // Construction of output times
    std::vector<Output::PairRepeatEachSteps> repeats_each_steps;

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

    bool const output_iteration_results =
        //! \ogs_file_param{prj__time_loop__output__output_iteration_results}
        config.getConfigParameter<bool>("output_iteration_results", false);

    return std::make_unique<Output>(output_directory, prefix, compress_output,
                                    data_mode, output_iteration_results,
                                    std::move(repeats_each_steps));
}

}  // namespace ProcessLib
