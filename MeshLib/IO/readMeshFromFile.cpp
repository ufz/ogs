/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "readMeshFromFile.h"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <boost/algorithm/string/erase.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"

#ifdef USE_PETSC
#include "MeshLib/IO/MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace
{
MeshLib::Mesh* readMeshFromFileSerial(const std::string& file_name,
                                      bool const compute_element_neighbors)
{
    if (BaseLib::hasFileExtension(".msh", file_name))
    {
        MeshLib::IO::Legacy::MeshIO meshIO;
        return meshIO.loadMeshFromFile(file_name);
    }

    if (BaseLib::hasFileExtension(".vtu", file_name))
    {
        return MeshLib::IO::VtuInterface::readVTUFile(
            file_name, compute_element_neighbors);
    }

    if (BaseLib::hasFileExtension(".vtk", file_name))
    {
        return MeshLib::IO::VtuInterface::readVTKFile(
            file_name, compute_element_neighbors);
    }

    ERR("readMeshFromFile(): Unknown mesh file format in file {:s}.",
        file_name);
    return nullptr;
}
}  // namespace

namespace MeshLib
{
namespace IO
{
MeshLib::Mesh* readMeshFromFile(const std::string& file_name,
                                bool const compute_element_neighbors)
{
#ifdef USE_PETSC
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 1)
    {
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        if (world_size > 1)
        {
            MeshLib::IO::NodePartitionedMeshReader read_pmesh(MPI_COMM_WORLD);
            const std::string file_name_base =
                BaseLib::dropFileExtension(file_name);
            return read_pmesh.read(file_name_base);
        }
        if (world_size == 1)
        {
            std::unique_ptr<Mesh> mesh{
                readMeshFromFileSerial(file_name, compute_element_neighbors)};
            return new MeshLib::NodePartitionedMesh(*mesh);
        }
        return nullptr;
    }
#endif
    return readMeshFromFileSerial(file_name, compute_element_neighbors);
}

}  // end namespace IO
}  // end namespace MeshLib
