/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
        explicit PairRepeatEachSteps(int c, int e) : repeat(c), each_steps(e) {}

        const int repeat;      //!< Apply \c each_steps \c repeat times.
        const int each_steps;  //!< Do output every \c each_steps timestep.
    };

public:
    Output(std::string output_directory, std::string prefix, std::string suffix,
           bool const compress_output, std::string const& data_mode,
           bool const output_nonlinear_iteration_results,
           std::vector<PairRepeatEachSteps> repeats_each_steps,
           std::vector<double>&& fixed_output_times,
           OutputDataSpecification&& output_data_specification,
           std::vector<std::string>&& mesh_names_for_output,
           std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

    //! TODO doc. Opens a PVD file for each process.
    void addProcess(ProcessLib::Process const& process, const int process_id);

    //! Writes output for the given \c process if it should be written in the
    //! given \c timestep.
    void doOutput(Process const& process, const int process_id,
                  const int timestep, const double t,
                  std::vector<GlobalVector*> const& x);

    //! Writes output for the given \c process if it has not been written yet.
    //! This method is intended for doing output after the last timestep in
    //! order to make sure that its results are written.
    void doOutputLastTimestep(Process const& process, const int process_id,
                              const int timestep, const double t,
                              std::vector<GlobalVector*> const& x);

    //! Writes output for the given \c process.
    //! This method will always write.
    //! It is intended to write output in error handling routines.
    void doOutputAlways(Process const& process, const int process_id,
                        const int timestep, const double t,
                        std::vector<GlobalVector*> const& x);

    //! Writes output for the given \c process.
    //! To be used for debug output after an iteration of the nonlinear solver.
    void doOutputNonlinearIteration(Process const& process,
                                    const int process_id, const int timestep,
                                    const double t,
                                    std::vector<GlobalVector*> const& x,
                                    const int iteration);

    std::vector<double> getFixedOutputTimes() {return fixed_output_times_;}

private:
    struct OutputFile;
    void outputBulkMesh(OutputFile const& output_file,
                        MeshLib::IO::PVDFile* const pvd_file,
                        MeshLib::Mesh const& mesh,
                        double const t) const;

private:
    std::string const output_directory_;
    std::string const output_file_prefix_;
    std::string const output_file_suffix_;

    //! Enables or disables zlib-compression of the output files.
    bool const output_file_compression_;

    //! Chooses vtk's data mode for output following the enumeration given in
    /// the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int const output_file_data_mode_;
    bool const output_nonlinear_iteration_results_;

    //! Describes after which timesteps to write output.
    std::vector<PairRepeatEachSteps> repeats_each_steps_;

    //! Given times that steps have to reach.
    std::vector<double> fixed_output_times_;

    std::multimap<Process const*, MeshLib::IO::PVDFile> process_to_pvd_file_;

    /**
     * Get the address of a PVDFile from corresponding to the given process.
     * @param process    Process.
     * @param process_id Process ID.
     * @param mesh_name_for_output mesh name for the output.
     * @return Address of a PVDFile.
     */
    MeshLib::IO::PVDFile* findPVDFile(Process const& process,
                                      const int process_id,
                                      std::string const& mesh_name_for_output);

    //! Determines if there should be output at the given \c timestep or \c t.
    bool shallDoOutput(int timestep, double const t);

    OutputDataSpecification const output_data_specification_;
    std::vector<std::string> mesh_names_for_output_;
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes_;
};

}  // namespace ProcessLib
