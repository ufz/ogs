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
 */

#include "VtuInterface.h"

#include <vtkNew.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#if defined(USE_PETSC) || defined(USE_MPI)
#include <vtkXMLPUnstructuredGridWriter.h>
#endif
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "BaseLib/Logging.h"

#include <boost/algorithm/string/erase.hpp>

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include "BaseLib/FileTools.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"

namespace MeshLib
{
namespace IO
{

VtuInterface::VtuInterface(const MeshLib::Mesh* mesh, int dataMode, bool compress) :
    mesh_(mesh), data_mode_(dataMode), use_compressor_(compress)
{
    if (data_mode_ == vtkXMLWriter::Ascii && compress)
    {
        WARN("Ascii data cannot be compressed, ignoring compression flag.");
    }
}

MeshLib::Mesh* VtuInterface::readVTUFile(std::string const &file_name)
{
    if (!BaseLib::IsFileExisting(file_name)) {
        ERR("File '{:s}' does not exist.", file_name);
        return nullptr;
    }

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();

    vtkUnstructuredGrid* vtkGrid = reader->GetOutput();
    if (vtkGrid->GetNumberOfPoints() == 0)
    {
        ERR("Mesh '{:s}' contains zero points.", file_name);
        return nullptr;
    }

    std::string const mesh_name (BaseLib::extractBaseNameWithoutExtension(file_name));
    return MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkGrid, mesh_name);
}

#ifdef USE_PETSC
std::string getVtuFileNameForPetscOutputWithoutExtension(
    std::string const& file_name)
{
    auto const file_name_extension = BaseLib::getFileExtension(file_name);
    if (file_name_extension != ".vtu")
    {
        OGS_FATAL("Expected a .vtu file for petsc output.");
    }

    auto const file_name_base = boost::erase_last_copy(file_name, ".vtu");
    auto basename = BaseLib::extractBaseName(file_name_base);

    // Replace dots to underscores since the pvtu writing function drops all
    // characters starting from a dot.
    std::replace(basename.begin(), basename.end(), '.', '_');

    // Restore the dirname if any.
    auto const dirname = BaseLib::extractPath(file_name_base);
    return BaseLib::joinPaths(dirname, basename);
}
#endif

bool VtuInterface::writeToFile(std::string const &file_name)
{
#ifdef USE_PETSC
    auto const vtu_file_name =
        getVtuFileNameForPetscOutputWithoutExtension(file_name);
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    return writeVTU<vtkXMLPUnstructuredGridWriter>(vtu_file_name + ".pvtu",
                                                   mpi_size, rank);
#else
    return writeVTU<vtkXMLUnstructuredGridWriter>(file_name);
#endif
}

int writeVtu(MeshLib::Mesh const& mesh, std::string const& file_name,
             int const data_mode)
{
    MeshLib::IO::VtuInterface writer(&mesh, data_mode);
    auto const result = writer.writeToFile(file_name);
    if (!result)
    {
        ERR("writeMeshToFile(): Could not write mesh to '{:s}'.", file_name);
        return -1;
    }
    return 0;
}


} // end namespace IO
} // end namespace MeshLib
