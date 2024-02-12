/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Output.h"

#include <cassert>
#include <exception>
#include <fstream>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/range/conversion.hpp>
#include <vector>

#include "AddProcessDataToMesh.h"
#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
void addBulkMeshPropertyToSubMesh(MeshLib::Mesh const& bulk_mesh,
                                  MeshLib::Mesh& sub_mesh,
                                  std::string const& property_name)
{
    if (bulk_mesh == sub_mesh)
    {
        return;
    }

    if (!bulk_mesh.getProperties().existsPropertyVector<double>(
            std::string_view(property_name)))
    {
        return;
    }
    auto const& bulk_mesh_property =
        *bulk_mesh.getProperties().getPropertyVector<double>(
            std::string_view(property_name));
    auto const mesh_item_type = bulk_mesh_property.getMeshItemType();

    std::string_view mesh_item_type_string =
        [mesh_item_type, &property_name]() -> std::string_view
    {
        switch (mesh_item_type)
        {
            case MeshLib::MeshItemType::Node:
                return MeshLib::getBulkIDString(MeshLib::MeshItemType::Node);
                break;
            case MeshLib::MeshItemType::Cell:
                return MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell);
                break;
            case MeshLib::MeshItemType::Edge:
                WARN(
                    "Property '{}' is assigned to edges. But mappings from the "
                    "bulk edges to submesh edges isn't implemented. Mesh item "
                    "type 'Edge' is not supported, only 'Node' and 'Cell' are "
                    "implemented at the moment.",
                    property_name);
                return "";
                break;
            case MeshLib::MeshItemType::Face:
                WARN(
                    "Property '{}' is assigned to faces. But mappings from the "
                    "bulk faces to submesh faces isn't implemented. Mesh item "
                    "type 'Face' is not supported, only 'Node' and 'Cell' are "
                    "implemented at the moment.",
                    property_name);
                return "";
                break;
            case MeshLib::MeshItemType::IntegrationPoint:
                WARN(
                    "Property '{}' is assigned to integration points. But "
                    "mappings from the bulk integration points to submesh "
                    "integration points isn't implemented. Mesh item type "
                    "'IntegrationPoint' is not supported, only 'Node' and "
                    "'Cell' are implemented at the moment.",
                    property_name);
                return "";
                break;
            default:
                OGS_FATAL("Unknown mesh item type '{}'.",
                          static_cast<int>(mesh_item_type));
                return "";
        }
    }();

    if (mesh_item_type_string.empty())
    {
        return;
    }
    if (!sub_mesh.getProperties().existsPropertyVector<std::size_t>(
            mesh_item_type_string, mesh_item_type, 1))
    {
        WARN(
            "The property {} is required for output on the mesh {}, but it "
            "doesn't exist.",
            mesh_item_type_string, sub_mesh.getName());
        return;
    }

    auto const& bulk_ids =
        *sub_mesh.getProperties().getPropertyVector<std::size_t>(
            mesh_item_type_string);

    auto const number_of_components =
        bulk_mesh_property.getNumberOfGlobalComponents();
    auto& sub_mesh_property = *MeshLib::getOrCreateMeshProperty<double>(
        sub_mesh, property_name, mesh_item_type, number_of_components);

    for (std::size_t sub_mesh_node_id = 0;
         sub_mesh_node_id < bulk_ids.getNumberOfTuples();
         ++sub_mesh_node_id)
    {
        auto const& bulk_id = bulk_ids[sub_mesh_node_id];
        for (std::remove_cv_t<decltype(number_of_components)> c = 0;
             c < number_of_components;
             ++c)
        {
            sub_mesh_property.getComponent(sub_mesh_node_id, c) =
                bulk_mesh_property.getComponent(bulk_id, c);
        }
    }
}

bool Output::isOutputProcess(const int process_id, const Process& process) const
{
    auto const is_last_process =
        process_id == static_cast<int>(_output_processes.size()) - 1;

    return process.isMonolithicSchemeUsed()
           // For the staggered scheme for the coupling, only the last
           // process, which gives the latest solution within a coupling
           // loop, is allowed to make output.
           || is_last_process;
}

Output::Output(std::unique_ptr<OutputFormat>&& output_format,
               bool const output_nonlinear_iteration_results,
               OutputDataSpecification&& output_data_specification,
               std::vector<std::string>&& mesh_names_for_output,
               std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
    : _output_format(std::move(output_format)),
      _output_nonlinear_iteration_results(output_nonlinear_iteration_results),
      _output_data_specification(std::move(output_data_specification)),
      _mesh_names_for_output(std::move(mesh_names_for_output)),
      _meshes(meshes)
{
}

void Output::addProcess(ProcessLib::Process const& process)
{
    _output_processes.push_back(process);
    if (_mesh_names_for_output.empty())
    {
        _mesh_names_for_output.push_back(process.getMesh().getName());
    }
}

void Output::doNotProjectFromBulkMeshToSubmeshes(
    std::string const& property_name,
    MeshLib::MeshItemType const mesh_item_type)
{
    _do_not_project_from_bulk_mesh_to_submeshes.emplace(property_name,
                                                        mesh_item_type);
}

void Output::outputMeshes(
    const int timestep, const double t, const int iteration,
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes)
    const
{
    if (_output_data_specification.output_variables.empty())
    {
        // special case: no output properties specified => output all properties
        for (auto const& mesh : meshes)
        {
            for (auto [name, property] : mesh.get().getProperties())
            {
                property->is_for_output = true;
            }
        }
    }
    else
    {
        for (auto const& mesh : meshes)
        {
            for (auto [name, property] : mesh.get().getProperties())
            {
                // special case: always output OGS_VERSION
                if (name == "OGS_VERSION")
                {
                    property->is_for_output = true;
                    continue;
                }

                property->is_for_output =
                    _output_data_specification.output_variables.contains(name);
            }
        }
    }
    _output_format->outputMeshes(timestep, t, iteration, meshes,
                                 _output_data_specification.output_variables);
}

MeshLib::Mesh const& Output::prepareSubmesh(
    std::string const& submesh_output_name, Process const& process,
    const int process_id, double const t,
    std::vector<GlobalVector*> const& xs) const
{
    auto& submesh = MeshLib::findMeshByName(_meshes.get(), submesh_output_name);

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
    auto const& property_names =
        bulk_mesh.getProperties().getPropertyVectorNames();

    // TODO Once all processes have been refactored to use the new residuum
    // assembly logic, the functionality of this lambda should be refactored.
    // Currently (Jan '23) there is a difference between the logic using this
    // lambda and doNotProjectFromBulkMeshToSubmeshes(): The latter is
    // considered regardless of submesh dimension.
    auto is_residuum_field = [](std::string const& name) -> bool
    {
        using namespace std::literals::string_view_literals;
        static constexpr std::string_view endings[] = {
            "FlowRate"sv, "heat_flux"sv, "MaterialForces"sv, "NodalForces"sv,
            "NodalForcesJump"sv};
        auto ends_with = [&](std::string_view const& ending)
        { return name.ends_with(ending); };
        return std::find_if(std::begin(endings), std::end(endings),
                            ends_with) != std::end(endings);
    };

    for (auto const& name : property_names)
    {
        if (_do_not_project_from_bulk_mesh_to_submeshes.contains(
                {name, MeshLib::MeshItemType::Node}))
        {
            // the projection is disabled regardless of mesh and submesh
            // dimension
            continue;
        }

        if (bulk_mesh.getDimension() == submesh.getDimension())
        {
            // omit the 'simple' transfer of the properties in the if condition
            // on submeshes with equal dimension to the bulk mesh
            // for those data extra assembly is required
            if (is_residuum_field(name))
            {
                continue;
            }
            addBulkMeshPropertyToSubMesh(bulk_mesh, submesh, name);
        }
        else
        {
            // For residuum based properties it is assumed that the lower
            // dimensional mesh is a boundary mesh!
            addBulkMeshPropertyToSubMesh(bulk_mesh, submesh, name);
        }
    }
    return submesh;
}

void Output::doOutputAlways(Process const& process,
                            const int process_id,
                            int const timestep,
                            const double t,
                            int const iteration,
                            std::vector<GlobalVector*> const& xs) const
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

    outputMeshes(timestep, t, iteration, std::move(output_meshes));

    INFO("[time] Output of timestep {:d} took {:g} s.", timestep,
         time_output.elapsed());
}

void Output::doOutput(Process const& process,
                      const int process_id,
                      int const timestep,
                      const double t,
                      int const iteration,
                      std::vector<GlobalVector*> const& xs) const
{
    if (isOutputStep(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, iteration, xs);
    }
#ifdef OGS_USE_INSITU
    // Note: last time step may be output twice: here and in
    // doOutputLastTimestep() which throws a warning.
    InSituLib::CoProcess(process.getMesh(), t, timestep, false,
                         _output_format->directory);
#endif
}

void Output::doOutputLastTimestep(Process const& process,
                                  const int process_id,
                                  int const timestep,
                                  const double t,
                                  int const iteration,
                                  std::vector<GlobalVector*> const& xs) const
{
    if (!isOutputStep(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, iteration, xs);
    }
#ifdef OGS_USE_INSITU
    InSituLib::CoProcess(process.getMesh(), t, timestep, true,
                         _output_format->directory);
#endif
}

void Output::doOutputNonlinearIteration(
    Process const& process, const int process_id, int const timestep,
    const double t, int const iteration,
    std::vector<GlobalVector*> const& xs) const
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

    std::string const output_file_name = _output_format->constructFilename(
        process.getMesh().getName(), timestep, t, iteration);

    std::string const output_file_path =
        BaseLib::joinPaths(_output_format->directory, output_file_name);

    DBUG("output iteration results to {:s}", output_file_path);

    if (dynamic_cast<OutputVTKFormat*>(_output_format.get()))
    {
        outputMeshVtk(
            output_file_path, process.getMesh(), _output_format->compression,
            dynamic_cast<OutputVTKFormat*>(_output_format.get())->data_mode);
    }
    else
    {
        DBUG("non-linear iterations can only written in Vtk/VTU format.");
    }
    INFO("[time] Output took {:g} s.", time_output.elapsed());
}

bool Output::isOutputStep(int const timestep, double const t) const
{
    return _output_data_specification.isOutputStep(timestep, t);
}

std::vector<std::string> Output::getFileNamesForOutput() const
{
    auto construct_filename = ranges::views::transform(
        [&](auto const& output_name)
        { return _output_format->constructFilename(output_name, 0, 0, 0); });

    return _mesh_names_for_output | construct_filename |
           ranges::to<std::vector>;
}

std::vector<double> calculateUniqueFixedTimesForAllOutputs(
    std::vector<Output> const& outputs)
{
    std::vector<double> fixed_times;
    for (auto const& output : outputs)
    {
        auto const& output_fixed_times = output.getFixedOutputTimes();
        fixed_times.insert(fixed_times.end(), output_fixed_times.begin(),
                           output_fixed_times.end());
    }
    BaseLib::makeVectorUnique(fixed_times);
    return fixed_times;
}

std::ostream& operator<<(std::ostream& os, Output const& output)
{
    os << "Output::_output_data_specification:\t"
       << output._output_data_specification;
    os << "Output::_output_format:\t" << *(output._output_format);
    return os;
}

}  // namespace ProcessLib
