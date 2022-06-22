/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "OutputFile.h"

#include <cassert>
#include <exception>
#include <fstream>
#include <vector>

namespace ProcessLib
{
/**
 * Get the address of a PVDFile corresponding to the given process.
 * @param process    Process.
 * @param process_id Process ID.
 * @param mesh_name_for_output mesh name for the output.
 * @param process_to_pvd_file a multimap that holds the PVD files associated
 * with each process.
 * @return Address of a PVDFile.
 */
MeshLib::IO::PVDFile& findPVDFile(
    Process const& process, const int process_id, std::string const& filename,
    std::multimap<Process const*, MeshLib::IO::PVDFile>& process_to_pvd_file)
{
    auto range = process_to_pvd_file.equal_range(&process);
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
                   MeshLib::IO::PVDFile& pvd_file, MeshLib::Mesh const& mesh,
                   double const t, int const timestep, int const iteration)
{
    auto const name =
        output_file.constructFilename(mesh.getName(), timestep, t, iteration);
    pvd_file.addVTUFile(name, t);

    auto const path = BaseLib::joinPaths(output_file.directory, name);
    outputMeshVtk(path, mesh, output_file.compression,
                  dynamic_cast<OutputVtkFormat const&>(output_file).data_mode);
}

std::string OutputVtkFormat::constructPVDName(
    std::string const& mesh_name) const
{
    return BaseLib::joinPaths(
        directory,
        BaseLib::constructFormattedFileName(prefix, mesh_name, 0, 0, 0) +
            ".pvd");
}

OutputFile::OutputFile(std::string const& directory, std::string const& prefix,
                       std::string const& suffix, bool const compression)
    : directory(directory),
      prefix(prefix),
      suffix(suffix),
      compression(compression)
{
}

std::string OutputVtkFormat::constructFilename(std::string mesh_name,
                                               int const timestep,
                                               double const t,
                                               int const iteration) const
{
    return BaseLib::constructFormattedFileName(prefix, mesh_name, timestep, t,
                                               iteration) +
           BaseLib::constructFormattedFileName(suffix, mesh_name, timestep, t,
                                               iteration) +
           ".vtu";
}

std::string OutputXDMFHDF5Format::constructFilename(std::string mesh_name,
                                                    int const timestep,
                                                    double const t,
                                                    int const iteration) const
{
    return BaseLib::constructFormattedFileName(prefix, mesh_name, timestep, t,
                                               iteration) +
           ".xdmf";
}

void OutputXDMFHDF5Format::outputMeshXdmf(
    std::set<std::string> const& output_variables,
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes,
    int const timestep, double const t, int const iteration)
{
    // \TODO (tm) Refactor to a dedicated VTKOutput and XdmfOutput
    // The XdmfOutput will create on construction the XdmfHdfWriter
    if (!mesh_xdmf_hdf_writer)
    {
        auto name = constructFilename(meshes[0].get().getName(), timestep, t,
                                      iteration);
        std::filesystem::path path(BaseLib::joinPaths(directory, name));
        mesh_xdmf_hdf_writer = std::make_unique<MeshLib::IO::XdmfHdfWriter>(
            std::move(meshes), path, timestep, t, output_variables, compression,
            n_files);
    }
    else
    {
        mesh_xdmf_hdf_writer->writeStep(t);
    };
}

void OutputVtkFormat::outputMeshes(
    const Process& process, const int process_id, const int timestep,
    const double t, const int iteration,
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes,
    [[maybe_unused]] std::set<std::string> const& output_variables)
{
    for (auto const& mesh : meshes)
    {
        auto const filename = constructPVDName(mesh.get().getName());
        auto& pvd_file =
            findPVDFile(process, process_id, filename, process_to_pvd_file);
        outputMeshVtk(*this, pvd_file, mesh, t, timestep, iteration);
    }
}

void OutputVtkFormat::addProcess(
    ProcessLib::Process const& process,
    std::vector<std::string> const& mesh_names_for_output)
{
    for (auto const& mesh_output_name : mesh_names_for_output)
    {
        auto const filename = constructPVDName(mesh_output_name);
        process_to_pvd_file.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(&process),
                                    std::forward_as_tuple(filename));
    }
}
}  // namespace ProcessLib
