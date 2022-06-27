/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Output.h"

#include <cassert>
#include <exception>
#include <fstream>
#include <vector>

#include "AddProcessDataToMesh.h"
#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
void addBulkMeshNodePropertyToSubMesh(MeshLib::Mesh const& bulk_mesh,
                                      MeshLib::Mesh& sub_mesh,
                                      std::string const& property_name)
{
    if (bulk_mesh == sub_mesh)
    {
        return;
    }
    if (!bulk_mesh.getProperties().existsPropertyVector<double>(
            property_name, MeshLib::MeshItemType::Node, 1))
    {
        return;
    }
    if (!sub_mesh.getProperties().existsPropertyVector<std::size_t>(
            "bulk_node_ids", MeshLib::MeshItemType::Node, 1))
    {
        return;
    }

    auto const& bulk_mesh_property =
        *bulk_mesh.getProperties().getPropertyVector<double>(property_name);
    auto const& bulk_node_ids =
        *sub_mesh.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids");

    auto& sub_mesh_property = *MeshLib::getOrCreateMeshProperty<double>(
        sub_mesh, property_name, MeshLib::MeshItemType::Node, 1);

    std::transform(std::begin(bulk_node_ids), std::end(bulk_node_ids),
                   std::begin(sub_mesh_property),
                   [&bulk_mesh_property](auto const id)
                   { return bulk_mesh_property[id]; });
}


bool Output::isOutputStep(int timestep, double const t) const
{
    auto const fixed_output_time = std::lower_bound(
        cbegin(_output_data_specification.fixed_output_times),
        cend(_output_data_specification.fixed_output_times), t);
    if ((fixed_output_time !=
         cend(_output_data_specification.fixed_output_times)) &&
        (std::abs(*fixed_output_time - t) <
         std::numeric_limits<double>::epsilon()))
    {
        return true;
    }

    int each_steps = 1;

    for (auto const& pair : _output_data_specification.repeats_each_steps)
    {
        each_steps = pair.each_steps;

        if (timestep > pair.repeat * each_steps)
        {
            timestep -= pair.repeat * each_steps;
        }
        else
        {
            break;
        }
    }

    return timestep % each_steps == 0;
}

bool Output::isOutputProcess(const int process_id, const Process& process) const
{
    if (!dynamic_cast<OutputVTKFormat*>(output_file.get()))
    {
        return process.isMonolithicSchemeUsed();
    }

    auto const number_of_pvd_files =
        dynamic_cast<OutputVTKFormat*>(output_file.get())
            ->process_to_pvd_file.size();
    auto const n_processes =
        static_cast<int>(number_of_pvd_files / _mesh_names_for_output.size());

    auto const is_last_process = process_id == n_processes - 1;

    return process.isMonolithicSchemeUsed()
           // For the staggered scheme for the coupling, only the last
           // process, which gives the latest solution within a coupling
           // loop, is allowed to make output.
           || is_last_process;
}

Output::Output(std::unique_ptr<OutputFile> output_file,
               bool const output_nonlinear_iteration_results,
               OutputDataSpecification&& output_data_specification,
               std::vector<std::string>&& mesh_names_for_output,
               std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
    : output_file(std::move(output_file)),
      _output_nonlinear_iteration_results(output_nonlinear_iteration_results),
      _output_data_specification(std::move(output_data_specification)),
      _mesh_names_for_output(mesh_names_for_output),
      _meshes(meshes)
{
}

void Output::addProcess(ProcessLib::Process const& process)
{
    if (_mesh_names_for_output.empty())
    {
        _mesh_names_for_output.push_back(process.getMesh().getName());
    }
    output_file->addProcess(process, _mesh_names_for_output);
}

void Output::outputMeshes(
    const Process& process, const int process_id, const int timestep,
    const double t, const int iteration,
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes)
{
    output_file->outputMeshes(process, process_id, timestep, t, iteration,
                              meshes,
                              _output_data_specification.output_variables);
}

MeshLib::Mesh const& Output::prepareSubmesh(
    std::string const& submesh_output_name, Process const& process,
    const int process_id, double const t,
    std::vector<GlobalVector*> const& xs) const
{
    auto& submesh = *BaseLib::findElementOrError(
        begin(_meshes), end(_meshes),
        [&submesh_output_name](auto const& m)
        { return m->getName() == submesh_output_name; },
        "Need mesh '" + submesh_output_name + "' for the output.");

    DBUG("Found {:d} nodes for output at mesh '{:s}'.",
         submesh.getNumberOfNodes(), submesh.getName());

    bool const output_secondary_variables = false;

    // TODO Under the assumption that xs.size() and submesh do not change during
    // the simulation, process output data should not be recreated every time,
    // but should rather be computed only once and stored for later reuse.
    auto const process_output_data =
        createProcessOutputData(process, xs.size(), submesh);

    addProcessDataToMesh(t, xs, process_id, process_output_data,
                         output_secondary_variables,
                         _output_data_specification);

    auto const& bulk_mesh = process.getMesh();
    auto const& node_property_names =
        bulk_mesh.getProperties().getPropertyVectorNames(
            MeshLib::MeshItemType::Node);
    for (auto const& name : node_property_names)
    {
        addBulkMeshNodePropertyToSubMesh(bulk_mesh, submesh, name);
    }
    return submesh;
}

void Output::doOutputAlways(Process const& process,
                            const int process_id,
                            int const timestep,
                            const double t,
                            int const iteration,
                            std::vector<GlobalVector*> const& xs)
{
    BaseLib::RunTime time_output;
    time_output.start();

    bool const output_secondary_variables = true;
    auto const process_output_data =
        createProcessOutputData(process, xs.size(), process.getMesh());

    // Need to add variables of process to mesh even if no output takes place.
    addProcessDataToMesh(t, xs, process_id, process_output_data,
                         output_secondary_variables,
                         _output_data_specification);

    if (!isOutputProcess(process_id, process))
    {
        return;
    }

    std::vector<std::reference_wrapper<const MeshLib::Mesh>> output_meshes;
    for (auto const& mesh_output_name : _mesh_names_for_output)
    {
        if (process.getMesh().getName() == mesh_output_name)
        {
            // process related output
            output_meshes.emplace_back(process.getMesh());
        }
        else
        {
            // mesh related output
            auto const& submesh =
                prepareSubmesh(mesh_output_name, process, process_id, t, xs);
            output_meshes.emplace_back(submesh);
        }
    }

    outputMeshes(process, process_id, timestep, t, iteration,
                 std::move(output_meshes));

    INFO("[time] Output of timestep {:d} took {:g} s.", timestep,
         time_output.elapsed());
}

void Output::doOutput(Process const& process,
                      const int process_id,
                      int const timestep,
                      const double t,
                      int const iteration,
                      std::vector<GlobalVector*> const& xs)
{
    if (isOutputStep(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, iteration, xs);
    }
#ifdef USE_INSITU
    // Note: last time step may be output twice: here and in
    // doOutputLastTimestep() which throws a warning.
    InSituLib::CoProcess(process.getMesh(), t, timestep, false,
                         output_file->directory);
#endif
}

void Output::doOutputLastTimestep(Process const& process,
                                  const int process_id,
                                  int const timestep,
                                  const double t,
                                  int const iteration,
                                  std::vector<GlobalVector*> const& xs)
{
    if (!isOutputStep(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, iteration, xs);
    }
#ifdef USE_INSITU
    InSituLib::CoProcess(process.getMesh(), t, timestep, true,
                         output_file->directory);
#endif
}

void Output::doOutputNonlinearIteration(Process const& process,
                                        const int process_id,
                                        int const timestep, const double t,
                                        int const iteration,
                                        std::vector<GlobalVector*> const& xs)
{
    if (!_output_nonlinear_iteration_results)
    {
        return;
    }

    BaseLib::RunTime time_output;
    time_output.start();

    bool const output_secondary_variable = true;
    auto const process_output_data =
        createProcessOutputData(process, xs.size(), process.getMesh());

    addProcessDataToMesh(t, xs, process_id, process_output_data,
                         output_secondary_variable, _output_data_specification);

    if (!isOutputProcess(process_id, process))
    {
        return;
    }

    std::string const output_file_name = output_file->constructFilename(
        process.getMesh().getName(), timestep, t, iteration);

    std::string const output_file_path =
        BaseLib::joinPaths(output_file->directory, output_file_name);

    DBUG("output iteration results to {:s}", output_file_path);

    if (dynamic_cast<OutputVTKFormat*>(output_file.get()))
    {
        outputMeshVtk(
            output_file_path, process.getMesh(), output_file->compression,
            dynamic_cast<OutputVTKFormat*>(output_file.get())->data_mode);
    }
    else
    {
        DBUG("non-linear iterations can only written in Vtk/VTU format.");
    }
    INFO("[time] Output took {:g} s.", time_output.elapsed());
}
}  // namespace ProcessLib
