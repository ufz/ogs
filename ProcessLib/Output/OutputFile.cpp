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
    if (output_file.type == ProcessLib::OutputType::vtk)
    {
        pvd_file.addVTUFile(name, t);
    }

    auto const path = BaseLib::joinPaths(output_file.directory, name);
    outputMeshVtk(path, mesh, output_file.compression, output_file.data_mode);
}

std::string OutputFile::constructPVDName(std::string const& mesh_name) const
{
    return BaseLib::joinPaths(
        directory,
        BaseLib::constructFormattedFileName(prefix, mesh_name, 0, 0, 0) +
            ".pvd");
}

OutputFile::OutputFile(std::string const& directory, OutputType const type,
                       std::string const& prefix, std::string const& suffix,
                       int const data_mode_, bool const compression_,
                       unsigned int const n_files)
    : directory(directory),
      prefix(prefix),
      suffix(suffix),
      type(type),
      data_mode(data_mode_),
      compression(compression_),
      n_files(n_files)
{
}

std::string OutputFile::constructFilename(std::string mesh_name,
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

void OutputFile::outputMeshXdmf(
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

}  // namespace ProcessLib
