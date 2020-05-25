/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Output.h"

#include <cassert>
#include <fstream>
#include <vector>

#include "BaseLib/Logging.h"

#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "ProcessLib/Process.h"

namespace
{
//! Converts a vtkXMLWriter's data mode string to an int. See
/// Output::output_file_data_mode_.
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
        "Unsupported vtk output file data mode '{:s}'. Expected Ascii, "
        "Binary, or Appended.",
        data_mode);
}

std::string constructPVDName(std::string const& output_directory,
                             std::string const& output_file_prefix,
                             int const process_id,
                             std::string const& mesh_name)
{
    return BaseLib::joinPaths(
        output_directory,
        BaseLib::constructFormattedFileName(output_file_prefix, mesh_name,
                                            process_id, 0, 0) +
            ".pvd");
}
}  // namespace

namespace ProcessLib
{
bool Output::shallDoOutput(int timestep, double const t)
{
    int each_steps = 1;

    for (auto const& pair : repeats_each_steps_)
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

    bool make_output = timestep % each_steps == 0;

    if (fixed_output_times_.empty())
    {
        return make_output;
    }

    const double specific_time = fixed_output_times_.back();
    const double zero_threshold = std::numeric_limits<double>::min();
    if (std::fabs(specific_time - t) < zero_threshold)
    {
        fixed_output_times_.pop_back();
        make_output = true;
    }

    return make_output;
}

Output::Output(std::string output_directory, std::string output_file_prefix,
               std::string output_file_suffix, bool const compress_output,
               std::string const& data_mode,
               bool const output_nonlinear_iteration_results,
               std::vector<PairRepeatEachSteps> repeats_each_steps,
               std::vector<double>&& fixed_output_times,
               OutputDataSpecification&& output_data_specification,
               std::vector<std::string>&& mesh_names_for_output,
               std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
    : output_directory_(std::move(output_directory)),
      output_file_prefix_(std::move(output_file_prefix)),
      output_file_suffix_(std::move(output_file_suffix)),
      output_file_compression_(compress_output),
      output_file_data_mode_(convertVtkDataMode(data_mode)),
      output_nonlinear_iteration_results_(output_nonlinear_iteration_results),
      repeats_each_steps_(std::move(repeats_each_steps)),
      fixed_output_times_(std::move(fixed_output_times)),
      output_data_specification_(std::move(output_data_specification)),
      mesh_names_for_output_(mesh_names_for_output),
      meshes_(meshes)
{
}

void Output::addProcess(ProcessLib::Process const& process,
                        const int process_id)
{
    if (mesh_names_for_output_.empty())
    {
        mesh_names_for_output_.push_back(process.getMesh().getName());
    }


    for (auto const& mesh_output_name : mesh_names_for_output_)
    {
        auto const filename =
            constructPVDName(output_directory_, output_file_prefix_,
                             process_id, mesh_output_name);
        process_to_pvd_file_.emplace(std::piecewise_construct,
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
    auto const filename =
        constructPVDName(output_directory_, output_file_prefix_, process_id,
                         mesh_name_for_output);
    auto range = process_to_pvd_file_.equal_range(&process);
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
            "The given process is not contained in the output"
            " configuration. Aborting.");
    }

    return pvd_file;
}

struct Output::OutputFile
{
    OutputFile(std::string const& directory, std::string const& prefix,
               std::string const& suffix, std::string const& mesh_name,
               int const process_id, int const timestep, double const t,
               int const data_mode_, bool const compression_)
        : name(BaseLib::constructFormattedFileName(prefix, mesh_name,
                                                   process_id, timestep, t) +
               BaseLib::constructFormattedFileName(suffix, mesh_name,
                                                   process_id, timestep, t) +
               ".vtu"),
          path(BaseLib::joinPaths(directory, name)),
          data_mode(data_mode_),
          compression(compression_)
    {
    }

    std::string const name;
    std::string const path;
    //! Chooses vtk's data mode for output following the enumeration given
    /// in the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation
    /// http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int const data_mode;

    //! Enables or disables zlib-compression of the output files.
    bool const compression;
};

void Output::outputBulkMesh(OutputFile const& output_file,
                            MeshLib::IO::PVDFile* const pvd_file,
                            MeshLib::Mesh const& mesh,
                            double const t) const
{
    DBUG("output to {:s}", output_file.path);

    pvd_file->addVTUFile(output_file.name, t);

    makeOutput(output_file.path, mesh, output_file.compression,
               output_file.data_mode);
}

void Output::doOutputAlways(Process const& process,
                            const int process_id,
                            const int timestep,
                            const double t,
                            std::vector<GlobalVector*> const& x)
{
    BaseLib::RunTime time_output;
    time_output.start();

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        dof_tables.push_back(&process.getDOFTable(i));
    }

    bool output_secondary_variable = true;
    // Need to add variables of process to vtu even no output takes place.
    addProcessDataToMesh(t, x, process_id, process.getMesh(), dof_tables,
                         process.getProcessVariables(process_id),
                         process.getSecondaryVariables(),
                         output_secondary_variable,
                         process.getIntegrationPointWriter(process.getMesh()),
                         output_data_specification_);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(process_to_pvd_file_.size() /
                                         mesh_names_for_output_.size()) -
                            1 ||
          process.isMonolithicSchemeUsed()))
    {
        return;
    }

    auto output_bulk_mesh = [&]() {
        outputBulkMesh(
            OutputFile(output_directory_, output_file_prefix_,
                       output_file_suffix_, process.getMesh().getName(),
                       process_id, timestep, t, output_file_data_mode_,
                       output_file_compression_),
            findPVDFile(process, process_id, process.getMesh().getName()),
            process.getMesh(), t);
    };

    for (auto const& mesh_output_name : mesh_names_for_output_)
    {
        if (process.getMesh().getName() == mesh_output_name)
        {
            output_bulk_mesh();
            continue;
        }
        auto& mesh = *BaseLib::findElementOrError(
            begin(meshes_), end(meshes_),
            [&mesh_output_name](auto const& m) {
                return m->getName() == mesh_output_name;
            },
            "Need mesh '" + mesh_output_name + "' for the output.");

        std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();
        DBUG("Found {:d} nodes for output at mesh '{:s}'.", nodes.size(),
             mesh.getName());

        MeshLib::MeshSubset mesh_subset(mesh, nodes);
        std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>
            mesh_dof_tables;
        mesh_dof_tables.reserve(x.size());
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            mesh_dof_tables.push_back(
                process.getDOFTable(i).deriveBoundaryConstrainedMap(
                    std::move(mesh_subset)));
        }
        std::vector<NumLib::LocalToGlobalIndexMap const*>
            mesh_dof_table_pointers;
        mesh_dof_table_pointers.reserve(mesh_dof_tables.size());
        transform(cbegin(mesh_dof_tables), cend(mesh_dof_tables),
                  back_inserter(mesh_dof_table_pointers),
                  [](std::unique_ptr<NumLib::LocalToGlobalIndexMap> const& p) {
                      return p.get();
                  });

        output_secondary_variable = false;
        addProcessDataToMesh(t, x, process_id, mesh, mesh_dof_table_pointers,
                             process.getProcessVariables(process_id),
                             process.getSecondaryVariables(),
                             output_secondary_variable,
                             process.getIntegrationPointWriter(mesh),
                             output_data_specification_);

        // TODO (TomFischer): add pvd support here. This can be done if the
        // output is mesh related instead of process related. This would also
        // allow for merging bulk mesh output and arbitrary mesh output.

        OutputFile const output_file{output_directory_,
                                     output_file_prefix_,
                                     output_file_suffix_,
                                     mesh.getName(),
                                     process_id,
                                     timestep,
                                     t,
                                     output_file_data_mode_,
                                     output_file_compression_};

        auto pvd_file = findPVDFile(process, process_id, mesh.getName());
        pvd_file->addVTUFile(output_file.name, t);

        DBUG("output to {:s}", output_file.path);

        makeOutput(output_file.path, mesh, output_file.compression,
                   output_file.data_mode);
    }
    INFO("[time] Output of timestep {:d} took {:g} s.", timestep,
         time_output.elapsed());
}

void Output::doOutput(Process const& process,
                      const int process_id,
                      const int timestep,
                      const double t,
                      std::vector<GlobalVector*> const& x)
{
    if (shallDoOutput(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, x);
    }
#ifdef USE_INSITU
    // Note: last time step may be output twice: here and in
    // doOutputLastTimestep() which throws a warning.
    InSituLib::CoProcess(process.getMesh(), t, timestep, false);
#endif
}

void Output::doOutputLastTimestep(Process const& process,
                                  const int process_id,
                                  const int timestep,
                                  const double t,
                                  std::vector<GlobalVector*> const& x)
{
    if (!shallDoOutput(timestep, t))
    {
        doOutputAlways(process, process_id, timestep, t, x);
    }
#ifdef USE_INSITU
    InSituLib::CoProcess(process.getMesh(), t, timestep, true);
#endif
}

void Output::doOutputNonlinearIteration(Process const& process,
                                        const int process_id,
                                        const int timestep, const double t,
                                        std::vector<GlobalVector*> const& x,
                                        const int iteration)
{
    if (!output_nonlinear_iteration_results_)
    {
        return;
    }

    BaseLib::RunTime time_output;
    time_output.start();

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        dof_tables.push_back(&process.getDOFTable(i));
    }

    bool const output_secondary_variable = true;
    addProcessDataToMesh(t, x, process_id, process.getMesh(), dof_tables,
                         process.getProcessVariables(process_id),
                         process.getSecondaryVariables(),
                         output_secondary_variable,
                         process.getIntegrationPointWriter(process.getMesh()),
                         output_data_specification_);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(process_to_pvd_file_.size()) - 1 ||
          process.isMonolithicSchemeUsed()))
    {
        return;
    }

    // Only check whether a process data is available for output.
    findPVDFile(process, process_id, process.getMesh().getName());

    std::string const output_file_name =
        BaseLib::constructFormattedFileName(output_file_prefix_,
                                            process.getMesh().getName(),
                                            process_id, timestep, t) +
        "_nliter_" + std::to_string(iteration) + ".vtu";
    std::string const output_file_path =
        BaseLib::joinPaths(output_directory_, output_file_name);

    DBUG("output iteration results to {:s}", output_file_path);

    INFO("[time] Output took {:g} s.", time_output.elapsed());

    makeOutput(output_file_path, process.getMesh(), output_file_compression_,
               output_file_data_mode_);
}
}  // namespace ProcessLib
