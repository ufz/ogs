/**
 * \file
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <filesystem>
#include <vtkXMLWriter.h>

namespace MeshLib {
class Mesh;

namespace IO
{

#ifdef USE_PETSC
std::string getVtuFileNameForPetscOutputWithoutExtension(
    std::string const& file_name);
#endif

/**
 * \brief Reads and writes VtkXMLUnstructuredGrid-files (vtu) to and from OGS data structures.
 * This class is currently not inherited from Writer because VTK will implement
 * writing to a string from 6.2 onwards.
 */
class VtuInterface final
{
public:
    /// Provide the mesh to write and set if compression should be used.
    explicit VtuInterface(const MeshLib::Mesh* mesh,
                          int dataMode = vtkXMLWriter::Appended,
                          bool compressed = false);

    /// Read an unstructured grid from a VTU file.
    /// \return The converted mesh or a nullptr if reading failed
    static MeshLib::Mesh* readVTUFile(std::string const &file_name);

    /// Read an unstructured grid from a legacy VTK file.
    /// Other data structures such as PolyData are ignored.
    /// \return The converted mesh or a nullptr if reading failed
    static MeshLib::Mesh* readVTKFile(std::string const& file_name);

    /// Writes the given mesh to file.
    /// \return True on success, false on error
    bool writeToFile(std::filesystem::path const& file_path);

    /// Writes the given mesh to vtu file.
    /// \param file_name      File name.
    /// \param num_partitions Number of partitions to be merged.
    /// \param rank the rank of the mpi process.
    /// \return True on success, false on error
    template <typename UnstructuredGridWriter>
    bool writeVTU(std::string const& file_name, const int num_partitions = 1,
                  const int rank = 1);

private:
    const MeshLib::Mesh* _mesh;
    int _data_mode;
    bool _use_compressor;
};

int writeVtu(MeshLib::Mesh const& mesh, std::string const& file_name,
             int const data_mode = vtkXMLWriter::Appended);

} // end namespace IO
} // end namespace MeshLib

#include "VtuInterface-impl.h"
