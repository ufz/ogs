/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
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

std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>
computeDofTablesForSubmesh(ProcessLib::Process const& process,
                           MeshLib::Mesh const& submesh,
                           std::size_t const num_processes)
{
    std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>> mesh_dof_tables;
    mesh_dof_tables.reserve(num_processes);

    for (std::size_t i = 0; i < num_processes; ++i)
    {
        mesh_dof_tables.push_back(
            process.getDOFTable(i).deriveBoundaryConstrainedMap(
                MeshLib::MeshSubset{submesh, submesh.getNodes()}));
    }

    return mesh_dof_tables;
}

std::vector<NumLib::LocalToGlobalIndexMap const*> toVectorOfConstPointers(
    std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table_pointers;

    dof_table_pointers.reserve(dof_tables.size());
    transform(cbegin(dof_tables), cend(dof_tables),
              back_inserter(dof_table_pointers),
              [](std::unique_ptr<NumLib::LocalToGlobalIndexMap> const& p)
              { return p.get(); });

    return dof_table_pointers;
}

void outputMeshVtk(ProcessLib::OutputFile const& output_file,
                   MeshLib::IO::PVDFile* const pvd_file,
                   MeshLib::Mesh const& mesh,
                   double const t)
{
    DBUG("output to {:s}", output_file.path);

    if (output_file.type == ProcessLib::OutputType::vtk)
    {
        pvd_file->addVTUFile(output_file.name, t);
    }

    ProcessLib::makeOutput(output_file.path, mesh, output_file.compression,
                           output_file.data_mode);
}
}  // namespace

namespace ProcessLib
{
bool Output::shallDoOutput(int timestep, double const t) const
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
bool Output::shallDoOutputStage2(const int process_id,
                                 const Process& process) const
{
    return process.isMonolithicSchemeUsed()
           // For the staggered scheme for the coupling, only the last process,
           // which gives the latest solution within a coupling loop, is allowed
           // to make output.
           || process_id == static_cast<int>(_process_to_pvd_file.size() /
                                             _mesh_names_for_output.size()) -
                                1;
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

// TODO return a reference.
MeshLib::IO::PVDFile* Output::findPVDFile(
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

    return pvd_file;
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

            auto* const pvd_file =
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

MeshLib::Mesh const& Output::prepareNonBulkMesh(
    std::string const& mesh_output_name, Process const& process,
    const int process_id, double const t,
    std::vector<GlobalVector*> const& xs) const
{
    auto& non_bulk_mesh = *BaseLib::findElementOrError(
        begin(_meshes), end(_meshes),
        [&mesh_output_name](auto const& m)
        { return m->getName() == mesh_output_name; },
        "Need mesh '" + mesh_output_name + "' for the output.");

    DBUG("Found {:d} nodes for output at mesh '{:s}'.",
         non_bulk_mesh.getNumberOfNodes(), non_bulk_mesh.getName());

    // TODO do not recreate everytime when doing output
    auto const mesh_dof_tables =
        computeDofTablesForSubmesh(process, non_bulk_mesh, xs.size());

    auto const mesh_dof_table_pointers =
        toVectorOfConstPointers(mesh_dof_tables);

    bool const output_secondary_variables = false;
    addProcessDataToMesh(
        t, xs, process_id, non_bulk_mesh, mesh_dof_table_pointers, process,
        output_secondary_variables, _output_data_specification);

    return non_bulk_mesh;
}

// TODO refactor this next.
void Output::doOutputAlways(Process const& process,
                            const int process_id,
                            int const timestep,
                            const double t,
                            int const iteration,
                            std::vector<GlobalVector*> const& xs)
{
    BaseLib::RunTime time_output;
    time_output.start();

    // Need to add variables of process to vtu even if no output takes place.
    bool const output_secondary_variables = true;
    addProcessDataToMesh(t, xs, process_id, process.getMesh(), process,
                         output_secondary_variables,
                         _output_data_specification);

    if (!shallDoOutputStage2(process_id, process))
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
            auto const& non_bulk_mesh = prepareNonBulkMesh(
                mesh_output_name, process, process_id, t, xs);
            output_meshes.push_back(non_bulk_mesh);
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
                      std::vector<GlobalVector*> const& x)
{
    if (shallDoOutput(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, iteration, x);
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
                                  std::vector<GlobalVector*> const& x)
{
    if (!shallDoOutput(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, iteration, x);
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
                                        std::vector<GlobalVector*> const& x)
{
    if (!_output_nonlinear_iteration_results)
    {
        return;
    }

    BaseLib::RunTime time_output;
    time_output.start();

    bool const output_secondary_variable = true;
    addProcessDataToMesh(t, x, process_id, process.getMesh(), process,
                         output_secondary_variable, _output_data_specification);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(_process_to_pvd_file.size()) - 1 ||
          process.isMonolithicSchemeUsed()))
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
    makeOutput(output_file_path, process.getMesh(), _output_file_compression,
               _output_file_data_mode);
    INFO("[time] Output took {:g} s.", time_output.elapsed());
}
}  // namespace ProcessLib
