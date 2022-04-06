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

struct OutputFile
{
    OutputFile(std::string const& directory, OutputType const type,
               std::string const& prefix, std::string const& suffix,
               std::string const& mesh_name, int const timestep, double const t,
               int const iteration, int const data_mode_,
               bool const compression_,
               std::set<std::string> const& outputnames,
               unsigned int const n_files)
        : name(constructFilename(type, prefix, suffix, mesh_name, timestep, t,
                                 iteration)),
          path(BaseLib::joinPaths(directory, name)),
          type(type),
          data_mode(data_mode_),
          compression(compression_),
          outputnames(outputnames),
          n_files(n_files)
    {
    }

    std::string const name;
    std::string const path;
    OutputType const type;

    //! Chooses vtk's data mode for output following the enumeration given
    /// in the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation
    /// http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int const data_mode;

    //! Enables or disables zlib-compression of the output files.
    bool const compression;
    std::set<std::string> outputnames;
    unsigned int n_files;

    static std::string constructFilename(OutputType const type,
                                         std::string prefix, std::string suffix,
                                         std::string mesh_name,
                                         int const timestep, double const t,
                                         int const iteration)
    {
        std::map<OutputType, std::string> filetype_to_extension = {
            {OutputType::vtk, "vtu"}, {OutputType::xdmf, "xdmf"}};

        try
        {
            std::string extension = filetype_to_extension.at(type);
            return BaseLib::constructFormattedFileName(prefix, mesh_name,
                                                       timestep, t, iteration) +
                   BaseLib::constructFormattedFileName(suffix, mesh_name,
                                                       timestep, t, iteration) +
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
};
}  // namespace ProcessLib

namespace
{
//! Converts a vtkXMLWriter's data mode string to an int. See
/// Output::_output_file_data_mode.
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

void outputMeshVtk(ProcessLib::OutputFile const& output_file,
                   MeshLib::IO::PVDFile& pvd_file,
                   MeshLib::Mesh const& mesh,
                   double const t)
{
    if (output_file.type == ProcessLib::OutputType::vtk)
    {
        pvd_file.addVTUFile(output_file.name, t);
    }

    outputMeshVtk(output_file.path, mesh, output_file.compression,
                  output_file.data_mode);
}
}  // namespace

namespace ProcessLib
{
bool Output::isOutputStep(int timestep, double const t) const
{
    auto const fixed_output_time = std::lower_bound(
        cbegin(_fixed_output_times), cend(_fixed_output_times), t);
    if ((fixed_output_time != cend(_fixed_output_times)) &&
        (std::abs(*fixed_output_time - t) <
         std::numeric_limits<double>::epsilon()))
    {
        return true;
    }

    int each_steps = 1;

    for (auto const& pair : _repeats_each_steps)
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

    if (timestep % each_steps == 0)
    {
        return true;
    }

    return false;
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

Output::Output(std::string directory, OutputType file_type,
               std::string file_prefix, std::string file_suffix,
               bool const compress_output, unsigned int const n_files,
               std::string const& data_mode,
               bool const output_nonlinear_iteration_results,
               std::vector<PairRepeatEachSteps> repeats_each_steps,
               std::vector<double>&& fixed_output_times,
               OutputDataSpecification&& output_data_specification,
               std::vector<std::string>&& mesh_names_for_output,
               std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
    : _output_directory(std::move(directory)),
      _output_file_type(file_type),
      _output_file_prefix(std::move(file_prefix)),
      _output_file_suffix(std::move(file_suffix)),
      _output_file_compression(compress_output),
      _n_files(n_files),
      _output_file_data_mode(convertVtkDataMode(data_mode)),
      _output_nonlinear_iteration_results(output_nonlinear_iteration_results),
      _repeats_each_steps(std::move(repeats_each_steps)),
      _fixed_output_times(std::move(fixed_output_times)),
      _output_data_specification(std::move(output_data_specification)),
      _mesh_names_for_output(mesh_names_for_output),
      _meshes(meshes)
{
    if (!std::is_sorted(cbegin(_fixed_output_times), cend(_fixed_output_times)))
    {
        OGS_FATAL(
            "Vector of fixed output time steps passed to the Output "
            "constructor must be sorted");
    }
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
            _output_directory, _output_file_prefix, mesh_output_name);
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
        _output_directory, _output_file_prefix, mesh_name_for_output);
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
    std::vector<std::reference_wrapper<const MeshLib::Mesh>>
        meshes,
    int const timestep,
    double const t)
{
    // \TODO (tm) Refactor to a dedicated VTKOutput and XdmfOutput
    // The XdmfOutput will create on construction the XdmfHdfWriter
    if (!_mesh_xdmf_hdf_writer)
    {
        std::filesystem::path path(output_file.path);
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
    if (_output_file_type == ProcessLib::OutputType::vtk)
    {
        for (auto const& mesh : meshes)
        {
            OutputFile const file(
                _output_directory, _output_file_type, _output_file_prefix,
                _output_file_suffix, mesh.get().getName(), timestep, t,
                iteration, _output_file_data_mode, _output_file_compression,
                _output_data_specification.output_variables, 1);

            auto& pvd_file =
                findPVDFile(process, process_id, mesh.get().getName());
            ::outputMeshVtk(file, pvd_file, mesh, t);
        }
    }
    else if (_output_file_type == ProcessLib::OutputType::xdmf)
    {
        std::string name = meshes[0].get().getName();
        OutputFile const file(
            _output_directory, _output_file_type, _output_file_prefix, "", name,
            timestep, t, iteration, _output_file_data_mode,
            _output_file_compression,
            _output_data_specification.output_variables, _n_files);

        outputMeshXdmf(std::move(file), std::move(meshes), timestep, t);
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
            output_meshes.push_back(process.getMesh());
        }
        else
        {
            // mesh related output
            auto const& submesh =
                prepareSubmesh(mesh_output_name, process, process_id, t, xs);
            output_meshes.push_back(submesh);
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
                         _output_directory);
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
                         _output_directory);
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

    std::string const output_file_name = OutputFile::constructFilename(
        _output_file_type, _output_file_prefix, _output_file_suffix,
        process.getMesh().getName(), timestep, t, iteration);

    std::string const output_file_path =
        BaseLib::joinPaths(_output_directory, output_file_name);

    DBUG("output iteration results to {:s}", output_file_path);
    outputMeshVtk(output_file_path, process.getMesh(), _output_file_compression,
                  _output_file_data_mode);
    INFO("[time] Output took {:g} s.", time_output.elapsed());
}
}  // namespace ProcessLib
