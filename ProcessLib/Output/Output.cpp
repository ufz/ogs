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

#ifndef _WIN32
#ifndef __APPLE__
#include <cfenv>
#endif  // __APPLE__
#endif  // _WIN32

#include "AddProcessDataToMesh.h"
#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
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

}  // namespace ProcessLib

namespace
{
//! Converts a vtkXMLWriter's data mode string to an int. See
/// Output::OutputFile::data_mode.
int convertVtkDataMode(std::string const& data_mode)
{
    if (data_mode == "Ascii")
    {
        return 0;
    }
    if (data_mode == "Binary")
    {
        return 1;
    }
    if (data_mode == "Appended")
    {
        return 2;
    }
    OGS_FATAL(
        "Unsupported vtk output file data mode '{:s}'. Expected Ascii, Binary, "
        "or Appended.",
        data_mode);
}

std::string constructPVDName(std::string const& output_directory,
                             std::string const& output_file_prefix,
                             std::string const& mesh_name)
{
    return BaseLib::joinPaths(output_directory,
                              BaseLib::constructFormattedFileName(
                                  output_file_prefix, mesh_name, 0, 0, 0) +
                                  ".pvd");
}

void outputMeshVtk(std::string const& file_name, MeshLib::Mesh const& mesh,
                   bool const compress_output, int const data_mode)
{
    DBUG("Writing output to '{:s}'.", file_name);

    // Store floating-point exception handling. Output of NaN's triggers
    // floating point exceptions. Because we are not debugging VTK (or other
    // libraries) at this point, the possibly set exceptions are temporary
    // disabled and restored before end of the function.
#ifndef _WIN32
#ifndef __APPLE__
    fenv_t fe_env;
    fegetenv(&fe_env);
    fesetenv(FE_DFL_ENV);  // Set default environment effectively disabling
                           // exceptions.
#endif                     //_WIN32
#endif                     //__APPLE__

    MeshLib::IO::VtuInterface vtu_interface(&mesh, data_mode, compress_output);
    vtu_interface.writeToFile(file_name);

    // Restore floating-point exception handling.
#ifndef _WIN32
#ifndef __APPLE__
    fesetenv(&fe_env);
#endif  //_WIN32
#endif  //__APPLE__
}

void outputMeshVtk(ProcessLib::Output::OutputFile const& output_file,
                   MeshLib::IO::PVDFile& pvd_file, MeshLib::Mesh const& mesh,
                   double const t, int const timestep, int const iteration)
{
    auto const name =
        output_file.constructFilename(mesh.getName(), timestep, t, iteration);
    if (output_file.type == ProcessLib::OutputType::vtk)
    {
        pvd_file.addVTUFile(name, t);
    }

    auto const path = BaseLib::joinPaths(output_file.directory, name);
    outputMeshVtk(path, mesh, output_file.compression, output_file.data_mode);
}
}  // namespace

namespace ProcessLib
{
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
    auto const n_processes = static_cast<int>(_process_to_pvd_file.size() /
                                              _mesh_names_for_output.size());

    auto const is_last_process = process_id == n_processes - 1;

    return process.isMonolithicSchemeUsed()
           // For the staggered scheme for the coupling, only the last process,
           // which gives the latest solution within a coupling loop, is allowed
           // to make output.
           || is_last_process;
}

Output::OutputFile::OutputFile(std::string const& directory,
                               OutputType const type, std::string const& prefix,
                               std::string const& suffix, int const data_mode_,
                               bool const compression_,
                               std::set<std::string> const& outputnames,
                               unsigned int const n_files)
    : directory(directory),
      prefix(prefix),
      suffix(suffix),
      type(type),
      data_mode(data_mode_),
      compression(compression_),
      outputnames(outputnames),
      n_files(n_files)
{
}

std::string Output::OutputFile::constructFilename(std::string mesh_name,
                                                  int const timestep,
                                                  double const t,
                                                  int const iteration) const
{
    std::map<OutputType, std::string> filetype_to_extension = {
        {OutputType::vtk, "vtu"}, {OutputType::xdmf, "xdmf"}};

    try
    {
        std::string extension = filetype_to_extension.at(type);
        return BaseLib::constructFormattedFileName(prefix, mesh_name, timestep,
                                                   t, iteration) +
               BaseLib::constructFormattedFileName(suffix, mesh_name, timestep,
                                                   t, iteration) +
               "." + extension;
    }
    catch (std::out_of_range&)
    {
        OGS_FATAL(
            "No supported file type provided. Read `{:s}' from <output><type> \
                in prj file. Supported: VTK, XDMF.",
            type);
    }
}

Output::Output(std::string directory, OutputType file_type,
               std::string file_prefix, std::string file_suffix,
               bool const compress_output, unsigned int const n_files,
               std::string const& data_mode,
               bool const output_nonlinear_iteration_results,
               OutputDataSpecification&& output_data_specification,
               std::vector<std::string>&& mesh_names_for_output,
               std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
    : output_file(directory, file_type, file_prefix, file_suffix,
                  convertVtkDataMode(data_mode), compress_output,
                  output_data_specification.output_variables, n_files),
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

    for (auto const& mesh_output_name : _mesh_names_for_output)
    {
        auto const filename = constructPVDName(
            output_file.directory, output_file.prefix, mesh_output_name);
        _process_to_pvd_file.emplace(std::piecewise_construct,
                                     std::forward_as_tuple(&process),
                                     std::forward_as_tuple(filename));
    }
}

MeshLib::IO::PVDFile& Output::findPVDFile(
    Process const& process,
    const int process_id,
    std::string const& mesh_name_for_output)
{
    auto const filename = constructPVDName(
        output_file.directory, output_file.prefix, mesh_name_for_output);
    auto range = _process_to_pvd_file.equal_range(&process);
    int counter = 0;
    MeshLib::IO::PVDFile* pvd_file = nullptr;
    for (auto spd_it = range.first; spd_it != range.second; ++spd_it)
    {
        if (spd_it->second.pvd_filename == filename)
        {
            if (counter == process_id)
            {
                pvd_file = &spd_it->second;
                break;
            }
            counter++;
        }
    }
    if (pvd_file == nullptr)
    {
        OGS_FATAL(
            "The given process is not contained in the output configuration. "
            "Aborting.");
    }

    return *pvd_file;
}

void Output::outputMeshXdmf(
    OutputFile const& output_file,
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes,
    int const timestep, double const t, int const iteration)
{
    // \TODO (tm) Refactor to a dedicated VTKOutput and XdmfOutput
    // The XdmfOutput will create on construction the XdmfHdfWriter
    if (!_mesh_xdmf_hdf_writer)
    {
        auto name = output_file.constructFilename(meshes[0].get().getName(),
                                                  timestep, t, iteration);
        std::filesystem::path path(
            BaseLib::joinPaths(output_file.directory, name));
        _mesh_xdmf_hdf_writer = std::make_unique<MeshLib::IO::XdmfHdfWriter>(
            std::move(meshes), path, timestep, t,
            _output_data_specification.output_variables,
            output_file.compression, output_file.n_files);
    }
    else
    {
        _mesh_xdmf_hdf_writer->writeStep(t);
    };
}

void Output::outputMeshes(
    const Process& process, const int process_id, const int timestep,
    const double t, const int iteration,
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes)
{
    if (output_file.type == ProcessLib::OutputType::vtk)
    {
        for (auto const& mesh : meshes)
        {
            OutputFile const file(
                output_file.directory, output_file.type, output_file.prefix,
                output_file.suffix, output_file.data_mode,
                output_file.compression,
                _output_data_specification.output_variables, 1);

            auto& pvd_file =
                findPVDFile(process, process_id, mesh.get().getName());
            ::outputMeshVtk(file, pvd_file, mesh, t, timestep, iteration);
        }
    }
    else if (output_file.type == ProcessLib::OutputType::xdmf)
    {
        std::string name = meshes[0].get().getName();
        OutputFile const file(
            output_file.directory, output_file.type, output_file.prefix, "",
            output_file.data_mode, output_file.compression,
            _output_data_specification.output_variables, output_file.n_files);

        outputMeshXdmf(std::move(file), std::move(meshes), timestep, t,
                       iteration);
    }
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
                         output_file.directory);
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
                         output_file.directory);
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

    // Only check whether a process data is available for output.
    findPVDFile(process, process_id, process.getMesh().getName());

    std::string const output_file_name = output_file.constructFilename(
        process.getMesh().getName(), timestep, t, iteration);

    std::string const output_file_path =
        BaseLib::joinPaths(output_file.directory, output_file_name);

    DBUG("output iteration results to {:s}", output_file_path);
    outputMeshVtk(output_file_path, process.getMesh(), output_file.compression,
                  output_file.data_mode);
    INFO("[time] Output took {:g} s.", time_output.elapsed());
}
}  // namespace ProcessLib
