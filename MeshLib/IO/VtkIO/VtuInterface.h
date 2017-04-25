/**
 * \file
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vtkXMLWriter.h>

namespace MeshLib {
class Mesh;

namespace IO
{

/**
 * \brief Reads and writes VtkXMLUnstructuredGrid-files (vtu) to and from OGS data structures.
 * This class is currently not inherited from Writer because VTK will implement
 * writing to a string from 6.2 onwards.
 */
class VtuInterface final
{
public:
    /// Provide the mesh to write and set if compression should be used.
    VtuInterface(const MeshLib::Mesh* mesh, int dataMode = vtkXMLWriter::Binary, bool compressed = false);

    /// Read an unstructured grid from a VTU file
    /// \return The converted mesh or a nullptr if reading failed
    static MeshLib::Mesh* readVTUFile(std::string const &file_name);

    /// Writes the given mesh to file.
    /// \return True on success, false on error
    bool writeToFile(std::string const &file_name);

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

} // end namespace IO
} // end namespace MeshLib

#include "VtuInterface-impl.h"
