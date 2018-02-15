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
    struct PairRepeatEachSteps
    {
        explicit PairRepeatEachSteps(unsigned c, unsigned e)
            : repeat(c), each_steps(e)
        {
        }

        const unsigned repeat;      //!< Apply \c each_steps \c repeat times.
        const unsigned each_steps;  //!< Do output every \c each_steps timestep.
    };

public:
    Output(std::string output_directory, std::string prefix,
           bool const compress_output, std::string const& data_mode,
           bool const output_nonlinear_iteration_results,
           std::vector<PairRepeatEachSteps> repeats_each_steps,
           const std::vector<double>&& specific_times);

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

    std::vector<double> getSpecificTimes() {return _specific_times;}

private:
    struct ProcessData
    {
        ProcessData(std::string const& filename) : pvd_file(filename) {}

        MeshLib::IO::PVDFile pvd_file;
    };

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

    //! Given times that steps have to reach.
    std::vector<double> _specific_times;

    std::multimap<Process const*, ProcessData> _process_to_process_data;

    /**
     * Get the address of a ProcessData from corresponding to the given process.
     * @param process    Process.
     * @param process_id Process ID.
     * @return Address of a ProcessData.
     */
    ProcessData* findProcessData(Process const& process, const int process_id);

    //! Determines if there should be output at the given \c timestep or \c t.
    bool shallDoOutput(unsigned timestep, double const t);
};



}  // namespace ProcessLib
