/**
 * \file VtuInterface-impl.h
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \author Wenqing Wang
 * \date   2015-10-19
 * \brief  Added parallel output of vtu pieces, and visualisation
 *         of the vtu picess via pvtu.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VtuInterface.h"

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <logog/include/logog.hpp>

#include "MeshLib/Vtk/VtkMappedMeshSource.h"

namespace MeshLib
{
namespace IO
{

template<typename UnstructuredGridWriter>
bool VtuInterface::writeVTU(std::string const &file_name, const int num_partitions)
{
    if(!_mesh)
    {
        ERR("VtuInterface::write(): No mesh specified.");
        return false;
    }

    vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
    vtkSource->SetMesh(_mesh);

    vtkSmartPointer<UnstructuredGridWriter> vtuWriter =
        vtkSmartPointer<UnstructuredGridWriter>::New();

    vtuWriter->SetInputConnection(vtkSource->GetOutputPort());

    if(_use_compressor)
        vtuWriter->SetCompressorTypeToZLib();
    else
        vtuWriter->SetCompressorTypeToNone();

    vtuWriter->SetDataMode(_data_mode);
    if (_data_mode == vtkXMLWriter::Appended)
        vtuWriter->SetEncodeAppendedData(1);
    if (_data_mode == vtkXMLWriter::Ascii)
    {
        // Mapped data structures for OGS to VTK mesh conversion are not fully
        // implemented and doing so is not trivial. Therefore for ascii output
        // the mapped unstructured grid is copied to a regular VTK grid.
        // See http://www.vtk.org/pipermail/vtkusers/2014-October/089400.html
        WARN(
            "Ascii data mode is currently not supported and the program may "
            "crash!");
        vtkSource->Update();
        vtkSmartPointer<vtkUnstructuredGrid> tempGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
        tempGrid->DeepCopy(vtkSource->GetOutput());
        vtuWriter->SetInputDataObject(tempGrid);
    }

    vtuWriter->SetFileName(file_name.c_str());
    if (num_partitions > 0)
        vtuWriter->SetNumberOfPieces(num_partitions);

    return (vtuWriter->Write() > 0);
}

} // end namespace IO
} // end namespace MeshLib
