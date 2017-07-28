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

#include <logog/include/logog.hpp>

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
    _mesh(mesh), _data_mode(dataMode), _use_compressor(compress)
{
    if(_data_mode == vtkXMLWriter::Ascii && compress)
        WARN("Ascii data cannot be compressed, ignoring compression flag.")
}

MeshLib::Mesh* VtuInterface::readVTUFile(std::string const &file_name)
{
    if (!BaseLib::IsFileExisting(file_name)) {
        ERR("File \"%s\" does not exist.", file_name.c_str());
        return nullptr;
    }

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();

    vtkUnstructuredGrid* vtkGrid = reader->GetOutput();
    if (vtkGrid->GetNumberOfPoints() == 0)
    {
        ERR("Mesh \"%s\" contains zero points.", file_name.c_str());
        return nullptr;
    }

    std::string const mesh_name (BaseLib::extractBaseNameWithoutExtension(file_name));
    return MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkGrid, mesh_name);
}

bool VtuInterface::writeToFile(std::string const &file_name)
{
#ifdef USE_PETSC
    auto const file_name_base = boost::erase_last_copy(file_name, ".vtu");

    auto const dirname = BaseLib::extractPath(file_name_base);
    auto basename = BaseLib::extractBaseName(file_name_base);

    // Since the pvtu writing function drops all letters from the letter of '.'.
    std::replace(basename.begin(), basename.end(), '.', '_');

    auto const vtu_file_name = BaseLib::joinPaths(dirname, basename);
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
} // end namespace IO
} // end namespace MeshLib
