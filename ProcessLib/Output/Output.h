/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <utility>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/IO/VtkIO/PVDFile.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"
#include "OutputDataSpecification.h"

namespace ProcessLib
{
class Process;

/// Writes output to the given \c file_name using the specified file format.
///
/// See Output::_output_file_data_mode documentation for the data_mode
/// parameter.
enum class OutputType : uint8_t
{
    vtk,
    xdmf
};

struct OutputFile
{
    OutputFile(std::string const& directory, OutputType const type,
               std::string const& prefix, std::string const& suffix,
               int const data_mode_, bool const compression_,
               std::set<std::string> const& outputnames,
               unsigned int const n_files);

    std::string const directory;
    std::string const prefix;
    std::string const suffix;
    OutputType const type;

    //! Chooses vtk's data mode for output following the enumeration given
    /// in the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation
    /// http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int const data_mode;

    //! Enables or disables zlib-compression of the output files.
    bool const compression;

    std::set<std::string> outputnames;

    std::unique_ptr<MeshLib::IO::XdmfHdfWriter> mesh_xdmf_hdf_writer;
    //! Specifies the number of hdf5 output files.
    unsigned int const n_files;

    std::string constructFilename(std::string mesh_name, int const timestep,
                                  double const t, int const iteration) const;
    void outputMeshXdmf(
        OutputDataSpecification const& output_data_specification,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes,
        int const timestep, double const t, int const iteration);
};

/*! Manages writing the solution of processes to disk.
 *
 * This class decides at which timesteps output is written
 * and initiates the writing process.
 */
class Output
{
public:
    Output(std::string directory, OutputType const file_type,
           std::string file_prefix, std::string file_suffix,
           bool const compress_output, unsigned int n_files,
           std::string const& data_mode,
           bool const output_nonlinear_iteration_results,
           OutputDataSpecification&& output_data_specification,
           std::vector<std::string>&& mesh_names_for_output,
           std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

    //! TODO doc. Opens a PVD file for each process.
    void addProcess(ProcessLib::Process const& process);

    //! Writes output for the given \c process if it should be written in the
    //! given \c timestep.
    void doOutput(Process const& process, const int process_id,
                  int const timestep, const double t, int const iteration,
                  std::vector<GlobalVector*> const& xs);

    //! Writes output for the given \c process if it has not been written yet.
    //! This method is intended for doing output after the last timestep in
    //! order to make sure that its results are written.
    void doOutputLastTimestep(Process const& process, const int process_id,
                              int const timestep, const double t,
                              int const iteration,
                              std::vector<GlobalVector*> const& xs);

    //! Writes output for the given \c process.
    //! This method will always write.
    //! It is intended to write output in error handling routines.
    void doOutputAlways(Process const& process, const int process_id,
                        int const timestep, const double t, int const iteration,
                        std::vector<GlobalVector*> const& xs);

    //! Writes output for the given \c process.
    //! To be used for debug output after an iteration of the nonlinear solver.
    void doOutputNonlinearIteration(Process const& process,
                                    const int process_id, int const timestep,
                                    const double t, const int iteration,
                                    std::vector<GlobalVector*> const& xs);

    std::vector<double> const& getFixedOutputTimes() const
    {
        return _output_data_specification.fixed_output_times;
    }

private:
    /**
     * Get the address of a PVDFile corresponding to the given process.
     * @param process    Process.
     * @param process_id Process ID.
     * @param mesh_name_for_output mesh name for the output.
     * @return Address of a PVDFile.
     */
    MeshLib::IO::PVDFile& findPVDFile(Process const& process,
                                      const int process_id,
                                      std::string const& mesh_name_for_output);

    //! Determines if there should be output at the given \c timestep or \c t.
    bool isOutputStep(int timestep, double const t) const;

    //! Determines if output should be written for the given process.
    //!
    //! With staggered coupling not every process writes output.
    bool isOutputProcess(int const process_id, Process const& process) const;

    void outputMeshes(
        Process const& process, const int process_id, int const timestep,
        const double t, int const iteration,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes);

    MeshLib::Mesh const& prepareSubmesh(
        std::string const& submesh_output_name, Process const& process,
        const int process_id, double const t,
        std::vector<GlobalVector*> const& xs) const;

    OutputFile output_file;

    bool const _output_nonlinear_iteration_results;

    //! Holds the PVD files associated with each process.
    //!
    //! Each \c process_id of each Process (in the current simulation) has a PVD
    //! file in this map for each element of #_mesh_names_for_output. I.e., the
    //! number of elements in this map is (roughly):
    //! <no. of processes> x <no. of process IDs per process> x <no. of meshes>.
    std::multimap<Process const*, MeshLib::IO::PVDFile> _process_to_pvd_file;

    OutputDataSpecification const _output_data_specification;
    std::vector<std::string> _mesh_names_for_output;
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& _meshes;
};

}  // namespace ProcessLib
