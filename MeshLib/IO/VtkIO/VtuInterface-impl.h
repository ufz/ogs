/**
 * \file
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
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VtuInterface.h"

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "BaseLib/Logging.h"

#include "MeshLib/Vtk/VtkMappedMeshSource.h"

namespace MeshLib
{
namespace IO
{
template <typename UnstructuredGridWriter>
bool VtuInterface::writeVTU(std::string const& file_name,
                            const int num_partitions, const int rank)
{
    if(!mesh_)
    {
        ERR("VtuInterface::write(): No mesh specified.");
        return false;
    }

    vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
    vtkSource->SetMesh(mesh_);

    vtkSmartPointer<UnstructuredGridWriter> vtuWriter =
        vtkSmartPointer<UnstructuredGridWriter>::New();

    vtkSource->Update();
    vtuWriter->SetInputData(vtkSource->GetOutput());

    if (use_compressor_)
    {
        vtuWriter->SetCompressorTypeToZLib();
    }
    else
    {
        vtuWriter->SetCompressorTypeToNone();
    }

    vtuWriter->SetDataMode(data_mode_);
    if (data_mode_ == vtkXMLWriter::Appended)
    {
        vtuWriter->SetEncodeAppendedData(1);
    }
    if (data_mode_ == vtkXMLWriter::Ascii)
    {
        vtkSource->Update();
        vtkSmartPointer<vtkUnstructuredGrid> tempGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
        tempGrid->DeepCopy(vtkSource->GetOutput());
        vtuWriter->SetInputDataObject(tempGrid);
    }

    vtuWriter->SetFileName(file_name.c_str());
#ifdef USE_PETSC
    vtuWriter->SetGhostLevel(1);
    vtuWriter->SetNumberOfPieces(num_partitions);
    vtuWriter->SetStartPiece(rank);
    vtuWriter->SetEndPiece(rank);
#else
    // avoid unused parameter warnings.
    (void)num_partitions;
    (void)rank;
#endif

    return (vtuWriter->Write() > 0);
}

} // end namespace IO
} // end namespace MeshLib
