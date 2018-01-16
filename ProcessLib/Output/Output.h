/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <utility>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/IO/VtkIO/PVDFile.h"
#include "ProcessOutput.h"

namespace ProcessLib
{
class Process;

/*! Manages writing the solution of processes to disk.
 *
 * This class decides at which timesteps output is written
 * and initiates the writing process.
 */
class Output
{
public:
    static std::unique_ptr<Output> newInstance(
        const BaseLib::ConfigTree& config, const std::string& output_directory);

    //! TODO doc. Opens a PVD file for each process.
    void addProcess(ProcessLib::Process const& process, const int process_id);

    //! Writes output for the given \c process if it should be written in the
    //! given \c timestep.
    void doOutput(Process const& process, const int process_id,
                  ProcessOutput const& process_output, unsigned timestep,
                  const double t, GlobalVector const& x);

    //! Writes output for the given \c process if it has not been written yet.
    //! This method is intended for doing output after the last timestep in
    //! order to make sure that its results are written.
    void doOutputLastTimestep(Process const& process, const int process_id,
                              ProcessOutput const& process_output,
                              unsigned timestep, const double t,
                              GlobalVector const& x);

    //! Writes output for the given \c process.
    //! This method will always write.
    //! It is intended to write output in error handling routines.
    void doOutputAlways(Process const& process, const int process_id,
                        ProcessOutput const& process_output, unsigned timestep,
                        const double t, GlobalVector const& x);

    //! Writes output for the given \c process.
    //! To be used for debug output after an iteration of the nonlinear solver.
    void doOutputNonlinearIteration(Process const& process,
                                    const int process_id,
                                    ProcessOutput const& process_output,
                                    const unsigned timestep, const double t,
                                    GlobalVector const& x,
                                    const unsigned iteration);

    struct PairRepeatEachSteps
    {
        explicit PairRepeatEachSteps(unsigned c, unsigned e)
            : repeat(c), each_steps(e)
        {
        }

        const unsigned repeat;      //!< Apply \c each_steps \c repeat times.
        const unsigned each_steps;  //!< Do output every \c each_steps timestep.
    };

private:
    struct SingleProcessData
    {
        SingleProcessData(std::string const& filename) : pvd_file(filename) {}

        MeshLib::IO::PVDFile pvd_file;
    };

    Output(std::string output_directory, std::string prefix,
           bool const compress_output, std::string const& data_mode,
           bool const output_nonlinear_iteration_results);

    std::string const _output_directory;
    std::string const _output_file_prefix;

    //! Enables or disables zlib-compression of the output files.
    bool const _output_file_compression;

    //! Chooses vtk's data mode for output following the enumeration given in
    /// the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int const _output_file_data_mode;
    bool const _output_nonlinear_iteration_results;

    //! Describes after which timesteps to write output.
    std::vector<PairRepeatEachSteps> _repeats_each_steps;

    std::multimap<Process const*, SingleProcessData> _single_process_data;

    /**
     * Get the address of a SingleProcessData from _single_process_data.
     * @param process    Process.
     * @param process_id Process ID.
     * @return Address of a SingleProcessData.
     */
    SingleProcessData* findSingleProcessData(Process const& process,
                                             const int process_id);
};

}  // namespace ProcessLib
