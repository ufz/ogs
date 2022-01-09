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
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "BaseLib/Logging.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"
#include "VtuInterface.h"

#ifdef USE_PETSC
#include <mpi.h>
#endif

class vtkXMLPUnstructuredGridWriter;

namespace MeshLib
{
namespace IO
{
template <typename UnstructuredGridWriter>
bool VtuInterface::writeVTU(std::string const& file_name,
                            [[maybe_unused]] const int num_partitions,
                            [[maybe_unused]] const int rank)
{
    if (!_mesh)
    {
        ERR("VtuInterface::write(): No mesh specified.");
        return false;
    }

    vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
    vtkSource->SetMesh(_mesh);

    vtkSmartPointer<UnstructuredGridWriter> vtuWriter =
        vtkSmartPointer<UnstructuredGridWriter>::New();

    vtkSource->Update();
    vtuWriter->SetInputData(vtkSource->GetOutput());

    if (_use_compressor)
    {
        vtuWriter->SetCompressorTypeToZLib();
    }
    else
    {
        vtuWriter->SetCompressorTypeToNone();
    }

    vtuWriter->SetDataMode(_data_mode);
    if (_data_mode == vtkXMLWriter::Appended)
    {
        vtuWriter->SetEncodeAppendedData(1);
    }
    if (_data_mode == vtkXMLWriter::Ascii)
    {
        vtkSource->Update();
        vtkSmartPointer<vtkUnstructuredGrid> tempGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
        tempGrid->DeepCopy(vtkSource->GetOutput());
        vtuWriter->SetInputDataObject(tempGrid);
    }

    vtuWriter->SetFileName(file_name.c_str());

    if constexpr (std::is_same_v<UnstructuredGridWriter,
                                 vtkXMLPUnstructuredGridWriter>)
    {
        vtuWriter->SetGhostLevel(1);
        vtuWriter->SetNumberOfPieces(num_partitions);
        vtuWriter->SetStartPiece(rank);
        vtuWriter->SetEndPiece(rank);
    }

#ifdef VTK_USE_64BIT_IDS
    vtuWriter->SetHeaderTypeToUInt64();
    // set SetIdTypeToInt64() as well?
#endif

    return (vtuWriter->Write() > 0);
}

}  // end namespace IO
}  // end namespace MeshLib
