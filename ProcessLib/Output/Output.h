/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iosfwd>
#include <map>
#include <utility>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "OutputDataSpecification.h"
#include "OutputFormat.h"

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
    Output(std::unique_ptr<OutputFormat>&& output_format,
           bool const output_nonlinear_iteration_results,
           OutputDataSpecification&& output_data_specification,
           std::vector<std::string>&& mesh_names_for_output,
           std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

    Output(Output const& other) = delete;
    Output(Output&& other) = default;
    Output& operator=(Output const& src) = delete;
    Output& operator=(Output&& src) = default;
    ~Output() = default;

    //! TODO doc. Opens a PVD file for each process.
    void addProcess(ProcessLib::Process const& process);

    //! Declare that the specified mesh property should not be projected from
    //! the bulk mesh to submeshes during submesh output.
    void doNotProjectFromBulkMeshToSubmeshes(
        std::string const& property_name,
        MeshLib::MeshItemType const mesh_item_type);

    //! Writes output for the given \c process if it should be written in
    //! the given \c timestep.
    void doOutput(Process const& process, const int process_id,
                  int const timestep, const double t, int const iteration,
                  std::vector<GlobalVector*> const& xs) const;

    //! Writes output for the given \c process if it has not been written yet.
    //! This method is intended for doing output after the last timestep in
    //! order to make sure that its results are written.
    void doOutputLastTimestep(Process const& process, const int process_id,
                              int const timestep, const double t,
                              int const iteration,
                              std::vector<GlobalVector*> const& xs) const;

    //! Writes output for the given \c process.
    //! This method will always write.
    //! It is intended to write output in error handling routines.
    void doOutputAlways(Process const& process, const int process_id,
                        int const timestep, const double t, int const iteration,
                        std::vector<GlobalVector*> const& xs) const;

    //! Writes output for the given \c process.
    //! To be used for debug output after an iteration of the nonlinear solver.
    void doOutputNonlinearIteration(Process const& process,
                                    const int process_id, int const timestep,
                                    const double t, const int iteration,
                                    std::vector<GlobalVector*> const& xs) const;

    //! Tells if output will be written at the specified timestep/time.
    bool isOutputStep(int const timestep, double const t) const;

    std::vector<double> const& getFixedOutputTimes() const
    {
        return _output_data_specification.fixed_output_times;
    }

    std::vector<std::string> getFileNamesForOutput() const;

    friend std::ostream& operator<<(std::ostream& os, Output const& output);

private:
    //! Determines if output should be written for the given process.
    //!
    //! With staggered coupling not every process writes output.
    bool isOutputProcess(int const process_id, Process const& process) const;

    void outputMeshes(
        int const timestep, const double t, int const iteration,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes)
        const;

    MeshLib::Mesh const& prepareSubmesh(
        std::string const& submesh_output_name, Process const& process,
        const int process_id, double const t,
        std::vector<GlobalVector*> const& xs) const;

    std::unique_ptr<OutputFormat> _output_format;

    bool _output_nonlinear_iteration_results;

    OutputDataSpecification _output_data_specification;
    std::vector<std::reference_wrapper<Process const>> _output_processes;
    std::vector<std::string> _mesh_names_for_output;
    // The reference wrapper enables the compiler to generate a move constructor
    // and move assignment operator
    // - another possibility would be to transform the
    // std::vector<std::unique_ptr<MeshLib::Mesh>> to
    // std::vector<MeshLib::Mesh*>
    // - this issue could also be solved by (globally) storing the meshes into a
    // std::vector<std::shared_ptr<MeshLib::Mesh>>
    std::reference_wrapper<std::vector<std::unique_ptr<MeshLib::Mesh>> const>
        _meshes;

    std::set<std::pair<std::string, MeshLib::MeshItemType>>
        _do_not_project_from_bulk_mesh_to_submeshes;
};

std::vector<double> calculateUniqueFixedTimesForAllOutputs(
    std::vector<Output> const& outputs);

}  // namespace ProcessLib
