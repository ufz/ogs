/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include<fstream>
#include "Output.h"

#include<cassert>
#include<vector>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"

namespace
{
//! Determines if there should be output at the given \c timestep.
template<typename CountsSteps>
bool shallDoOutput(unsigned timestep, CountsSteps const& repeats_each_steps)
{
    unsigned each_steps = 1;

    for (auto const& pair : repeats_each_steps)
    {
        each_steps = pair.each_steps;

        if (timestep > pair.repeat * each_steps)
        {
            timestep -= pair.repeat * each_steps;
        } else {
            break;
        }
    }

    return timestep % each_steps == 0;
}
}

namespace ProcessLib
{

template<typename GlobalSetup>
std::unique_ptr<Output<GlobalSetup>> Output<GlobalSetup>::
newInstance(const BaseLib::ConfigTree &config, std::string const& output_directory)
{
    std::unique_ptr<Output> out{ new Output{
        BaseLib::joinPaths(output_directory,
            //! \ogs_file_param{prj__output__prefix}
            config.getParameter<std::string>("prefix"))}};

    //! \ogs_file_param{prj__output__timesteps}
    if (auto const timesteps = config.getSubtreeOptional("timesteps"))
    {
        //! \ogs_file_param{prj__output__timesteps__pair}
        for (auto pair : timesteps->getSubtreeList("pair"))
        {
            //! \ogs_file_param{prj__output__timesteps__pair__repeat}
            auto repeat     = pair.getParameter<unsigned>("repeat");
            //! \ogs_file_param{prj__output__timesteps__pair__each_steps}
            auto each_steps = pair.getParameter<unsigned>("each_steps");

            assert(repeat != 0 && each_steps != 0);
            out->_repeats_each_steps.emplace_back(repeat, each_steps);
        }

        if (out->_repeats_each_steps.empty()) {
            ERR("You have not given any pair (<repeat/>, <each_steps/>) that defines"
                    " at which timesteps output shall be written. Aborting.");
            std::abort();
        }
    }
    else
    {
        out->_repeats_each_steps.emplace_back(1, 1);
    }

    return out;
}

template<typename GlobalSetup>
void Output<GlobalSetup>::
initialize(typename Output<GlobalSetup>::ProcessIter first,
     typename Output<GlobalSetup>::ProcessIter last)
{
    for (unsigned pcs_idx = 0; first != last; ++first, ++pcs_idx)
    {
        auto const filename = _output_file_prefix
                              + "_pcs_" + std::to_string(pcs_idx)
                              + ".pvd";
        _single_process_data.emplace(std::piecewise_construct,
                std::forward_as_tuple(&**first),
                std::forward_as_tuple(pcs_idx, filename));
    }
}

template<typename GlobalSetup>
void Output<GlobalSetup>::
doOutputAlways(Process<GlobalSetup> const& process,
               unsigned timestep,
               const double t,
               typename GlobalSetup::VectorType const& x)
{
    auto spd_it = _single_process_data.find(&process);
    if (spd_it == _single_process_data.end()) {
        ERR("The given process is not contained in the output configuration."
            " Aborting.");
        std::abort();
    }
    auto& spd = spd_it->second;

    std::string const output_file_name =
            _output_file_prefix + "_pcs_" + std::to_string(spd.process_index)
            + "_ts_" + std::to_string(timestep)
            + "_t_"  + std::to_string(t)
            + ".vtu";
    DBUG("output to %s", output_file_name.c_str());
    process.output(output_file_name, timestep, x);
    spd.pvd_file.addVTUFile(output_file_name, t);
}

template<typename GlobalSetup>
void Output<GlobalSetup>::
doOutput(Process<GlobalSetup> const& process,
         unsigned timestep,
         const double t,
         typename GlobalSetup::VectorType const& x)
{
    if (shallDoOutput(timestep, _repeats_each_steps))
        doOutputAlways(process, timestep, t, x);
}

template<typename GlobalSetup>
void Output<GlobalSetup>::
doOutputLastTimestep(Process<GlobalSetup> const& process,
                     unsigned timestep,
                     const double t,
                     typename GlobalSetup::VectorType const& x)
{
    if (!shallDoOutput(timestep, _repeats_each_steps))
        doOutputAlways(process, timestep, t, x);
}


// explicit instantiation
template class Output<GlobalSetupType>;

}
