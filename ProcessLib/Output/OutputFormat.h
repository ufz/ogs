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
#include <vector>

#include "MeshLib/IO/VtkIO/PVDFile.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"

namespace ProcessLib
{
class Process;

struct OutputFormat
{
    OutputFormat(std::string const& directory, std::string prefix,
                 std::string suffix, bool const compression);
    virtual ~OutputFormat() = default;

    OutputFormat(OutputFormat const& other) = delete;
    OutputFormat(OutputFormat&& other) = default;
    OutputFormat& operator=(OutputFormat& other) = delete;
    OutputFormat& operator=(OutputFormat&& other) = default;

    std::string directory;
    std::string prefix;
    std::string suffix;
    //! Enables or disables zlib-compression of the output files.
    bool compression;

    virtual void outputMeshes(
        const int timestep, const double t, const int iteration,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
        std::set<std::string> const& output_variables) const = 0;
    virtual std::string constructFilename(std::string const& mesh_name,
                                          int const timestep,
                                          double const t,
                                          int const iteration) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, OutputFormat const& of)
{
    os << "OutputFormat::directory:" << of.directory << std::endl;
    os << "OutputFormat::prefix:" << of.prefix << std::endl;
    os << "OutputFormat::suffix:" << of.suffix << std::endl;
    os << "OutputFormat::compression:" << of.compression << std::endl;
    return os;
}

struct OutputVTKFormat final : public OutputFormat
{
    OutputVTKFormat(std::string const& directory, std::string prefix,
                    std::string suffix, bool const compression,
                    int const data_mode)
        : OutputFormat(directory, std::move(prefix), std::move(suffix),
                       compression),
          data_mode(data_mode)
    {
    }

    void outputMeshes(
        const int timestep, const double t, const int iteration,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
        std::set<std::string> const& output_variables) const override;

    //! Chooses vtk's data mode for output following the enumeration given
    /// in the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation
    /// http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int data_mode;

    //! Holds the PVD files associated with the meshes.
    mutable std::map<std::string, MeshLib::IO::PVDFile> mesh_name_to_pvd_file;

    MeshLib::IO::PVDFile& findOrCreatePVDFile(
        std::string const& mesh_name) const;

    std::string constructFilename(std::string const& mesh_name,
                                  int const timestep, double const t,
                                  int const iteration) const override;

    std::string constructPVDName(std::string const& mesh_name) const;
};

struct OutputXDMFHDF5Format final : public OutputFormat
{
    OutputXDMFHDF5Format(std::string const& directory, std::string prefix,
                         std::string suffix, bool const compression,
                         unsigned int const n_files,
                         unsigned int const chunk_size_bytes)
        : OutputFormat(directory, std::move(prefix), std::move(suffix),
                       compression),
          n_files(n_files),
          chunk_size_bytes(chunk_size_bytes)
    {
    }

    void outputMeshes(
        const int timestep, const double t, const int iteration,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
        std::set<std::string> const& output_variables) const override
    {
        outputMeshXdmf(output_variables, meshes, timestep, t, iteration);
    }

    std::string constructFilename(std::string const& mesh_name,
                                  int const timestep, double const t,
                                  int const iteration) const override;

    mutable std::unique_ptr<MeshLib::IO::XdmfHdfWriter> mesh_xdmf_hdf_writer;
    //! Specifies the number of hdf5 output files.
    unsigned int n_files;
    //! Specifies the chunks size in bytes per hdf5 output file.
    unsigned int const chunk_size_bytes;

    void outputMeshXdmf(
        std::set<std::string> const& output_variables,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
        int const timestep, double const t, int const iteration) const;
};

void outputMeshVtk(std::string const& file_name, MeshLib::Mesh const& mesh,
                   bool const compress_output, int const data_mode);
}  // namespace ProcessLib
