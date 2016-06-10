/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_OUTPUT_H
#define PROCESSLIB_OUTPUT_H

#include "BaseLib/ConfigTree.h"
#include "MeshLib/IO/VtkIO/PVDFile.h"
#include "Process.h"

namespace ProcessLib
{

/*! Manages writing the solution of processes to disk.
 *
 * This class decides at which timesteps output is written
 * and initiates the writing process.
 */
template<typename GlobalSetup>
class Output
{
public:
    static std::unique_ptr<Output>
    newInstance(const BaseLib::ConfigTree& config,
                const std::string& output_directory);

    using ProcessIter = std::vector<std::unique_ptr<ProcessLib::Process<GlobalSetupType>>>
                        ::const_iterator;

    //! Opens a PVD file for each process.
    void initialize(ProcessIter first, ProcessIter last);

    //! Writes output for the given \c process if it should be written in the
    //! given \c timestep.
    void doOutput(
            Process<GlobalSetup> const& process, unsigned timestep,
            const double t,
            typename GlobalSetup::VectorType const& x);

    //! Writes output for the given \c process if it has not been written yet.
    //! This method is intended for doing output after the last timestep in order
    //! to make sure that its results are written.
    void doOutputLastTimestep(
            Process<GlobalSetup> const& process, unsigned timestep,
            const double t,
            typename GlobalSetup::VectorType const& x);

    //! Writes output for the given \c process.
    //! This method will always write.
    //! It is intended to write output in error handling routines.
    void doOutputAlways(
            Process<GlobalSetup> const& process, unsigned timestep,
            const double t,
            typename GlobalSetup::VectorType const& x);

    struct PairRepeatEachSteps
    {
        explicit PairRepeatEachSteps(unsigned c, unsigned e)
            : repeat(c), each_steps(e)
        {}

        const unsigned repeat;     //!< Apply \c each_steps \c repeat times.
        const unsigned each_steps; //!< Do output every \c each_steps timestep.
    };
private:
    struct SingleProcessData
    {
        SingleProcessData(unsigned process_index_,
                          std::string const& filename)
            : process_index(process_index_)
            , pvd_file(filename)
        {}

        const unsigned process_index;
        MeshLib::IO::PVDFile pvd_file;
    };

    Output(std::string const& prefix)
        : _output_file_prefix(prefix)
    {}

    std::string _output_file_prefix;

    //! Describes after which timesteps to write output.
    std::vector<PairRepeatEachSteps> _repeats_each_steps;

    std::map<Process<GlobalSetup> const*, SingleProcessData> _single_process_data;
};

}

#endif // PROCESSLIB_OUTPUT_H
