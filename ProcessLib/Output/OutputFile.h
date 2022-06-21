/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MeshLib/IO/VtkIO/PVDFile.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"

namespace ProcessLib
{
enum class OutputType : uint8_t
{
    vtk,
    xdmf
};

struct OutputFile
{
    OutputFile(std::string const& directory, OutputType const type,
               std::string const& prefix, std::string const& suffix,
               int const data_mode_, bool const compression_,
               unsigned int const n_files);

    std::string directory;
    std::string prefix;
    std::string suffix;
    OutputType const type;

    //! Chooses vtk's data mode for output following the enumeration given
    /// in the vtkXMLWriter: {Ascii, Binary, Appended}.  See vtkXMLWriter
    /// documentation
    /// http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html
    int const data_mode;

    //! Enables or disables zlib-compression of the output files.
    bool const compression;

    std::unique_ptr<MeshLib::IO::XdmfHdfWriter> mesh_xdmf_hdf_writer;
    //! Specifies the number of hdf5 output files.
    unsigned int const n_files;

    std::string constructFilename(std::string mesh_name, int const timestep,
                                  double const t, int const iteration) const;

    void outputMeshXdmf(
        std::set<std::string> const& output_variables,
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes,
        int const timestep, double const t, int const iteration);

    std::string constructPVDName(std::string const& mesh_name) const;
};

void outputMeshVtk(std::string const& file_name, MeshLib::Mesh const& mesh,
                   bool const compress_output, int const data_mode);
}  // namespace ProcessLib
