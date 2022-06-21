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
    OutputDataSpecification const& output_data_specification,
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
            std::move(meshes), path, timestep, t,
            output_data_specification.output_variables, compression, n_files);
    }
    else
    {
        mesh_xdmf_hdf_writer->writeStep(t);
    };
}

}  // namespace ProcessLib
