/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Output.h"

#include <cassert>
#include <fstream>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "Applications/InSituLib/Adaptor.h"

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

std::unique_ptr<Output> Output::
newInstance(const BaseLib::ConfigTree &config, std::string const& output_directory)
{
    auto const output_iteration_results =
        //! \ogs_file_param{prj__time_loop__output__output_iteration_results}
        config.getConfigParameterOptional<bool>("output_iteration_results");

    std::unique_ptr<Output> out{new Output{
        output_directory,
        //! \ogs_file_param{prj__time_loop__output__prefix}
        config.getConfigParameter<std::string>("prefix"),
        //! \ogs_file_param{prj__time_loop__output__compress_output}
        config.getConfigParameter("compress_output", true),
        output_iteration_results ? *output_iteration_results : false}};

    //! \ogs_file_param{prj__time_loop__output__timesteps}
    if (auto const timesteps = config.getConfigSubtreeOptional("timesteps"))
    {
        //! \ogs_file_param{prj__time_loop__output__timesteps__pair}
        for (auto pair : timesteps->getConfigSubtreeList("pair"))
        {
            //! \ogs_file_param{prj__time_loop__output__timesteps__pair__repeat}
            auto repeat     = pair.getConfigParameter<unsigned>("repeat");
            //! \ogs_file_param{prj__time_loop__output__timesteps__pair__each_steps}
            auto each_steps = pair.getConfigParameter<unsigned>("each_steps");

            assert(repeat != 0 && each_steps != 0);
            out->_repeats_each_steps.emplace_back(repeat, each_steps);
        }

        if (out->_repeats_each_steps.empty()) {
            OGS_FATAL("You have not given any pair (<repeat/>, <each_steps/>) that defines"
                    " at which timesteps output shall be written. Aborting.");
        }
    }
    else
    {
        out->_repeats_each_steps.emplace_back(1, 1);
    }

    return out;
}

void Output::addProcess(ProcessLib::Process const& process, const unsigned pcs_idx)
{
    auto const filename =
        BaseLib::joinPaths(_output_directory, _output_file_prefix + "_pcs_" + std::to_string(pcs_idx) + ".pvd");
    _single_process_data.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(&process),
                                 std::forward_as_tuple(pcs_idx, filename));
}

void Output::doOutputAlways(Process const& process,
                            ProcessOutput const& process_output,
                            unsigned timestep,
                            const double t,
                            GlobalVector const& x)
{
    BaseLib::RunTime time_output;
    time_output.start();

    auto spd_it = _single_process_data.find(&process);
    if (spd_it == _single_process_data.end()) {
        OGS_FATAL("The given process is not contained in the output configuration."
            " Aborting.");
    }
    auto& spd = spd_it->second;

    std::string const output_file_name =
            _output_file_prefix + "_pcs_" + std::to_string(spd.process_index)
            + "_ts_" + std::to_string(timestep)
            + "_t_"  + std::to_string(t)
            + ".vtu";
    std::string const output_file_path = BaseLib::joinPaths(_output_directory, output_file_name);
    DBUG("output to %s", output_file_path.c_str());
    doProcessOutput(output_file_path, _output_file_compression, t, x,
                    process.getMesh(), process.getDOFTable(),
                    process.getProcessVariables(),
                    process.getSecondaryVariables(), process_output);
    spd.pvd_file.addVTUFile(output_file_name, t);

    INFO("[time] Output of timestep %d took %g s.", timestep,
         time_output.elapsed());
}

void Output::doOutput(Process const& process,
                      ProcessOutput const& process_output,
                      unsigned timestep,
                      const double t,
                      GlobalVector const& x)
{
    if (shallDoOutput(timestep, _repeats_each_steps))
        doOutputAlways(process, process_output, timestep, t, x);
#ifdef USE_INSITU
    // Note: last time step may be output twice: here and in
    // doOutputLastTimestep() which throws a warning.
    InSituLib::CoProcess(process.getMesh(), t, timestep, false);
#endif
}

void Output::doOutputLastTimestep(Process const& process,
                                  ProcessOutput const& process_output,
                                  unsigned timestep,
                                  const double t,
                                  GlobalVector const& x)
{
    if (!shallDoOutput(timestep, _repeats_each_steps))
        doOutputAlways(process, process_output, timestep, t, x);
#ifdef USE_INSITU
    InSituLib::CoProcess(process.getMesh(), t, timestep, true);
#endif
}

void Output::doOutputNonlinearIteration(Process const& process,
                                        ProcessOutput const& process_output,
                                        const unsigned timestep, const double t,
                                        GlobalVector const& x,
                                        const unsigned iteration) const
{
    if (!_output_nonlinear_iteration_results) return;

    BaseLib::RunTime time_output;
    time_output.start();

    auto spd_it = _single_process_data.find(&process);
    if (spd_it == _single_process_data.end()) {
        OGS_FATAL("The given process is not contained in the output configuration."
            " Aborting.");
    }
    auto& spd = spd_it->second;

    std::string const output_file_name =
            _output_file_prefix + "_pcs_" + std::to_string(spd.process_index)
            + "_ts_" + std::to_string(timestep)
            + "_t_"  + std::to_string(t)
            + "_nliter_" + std::to_string(iteration)
            + ".vtu";
    std::string const output_file_path = BaseLib::joinPaths(_output_directory, output_file_name);
    DBUG("output iteration results to %s", output_file_path.c_str());
    doProcessOutput(output_file_path, _output_file_compression, t, x,
                    process.getMesh(), process.getDOFTable(),
                    process.getProcessVariables(),
                    process.getSecondaryVariables(), process_output);

    INFO("[time] Output took %g s.", time_output.elapsed());
}

}
