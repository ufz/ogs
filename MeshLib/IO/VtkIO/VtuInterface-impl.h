/**
 * \file
 * \author Lars Bilke
 * \author Wenqing Wang
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "BaseLib/Logging.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"
#include "VtuInterface.h"

#ifdef USE_PETSC
#include <mpi.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
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

#ifdef USE_PETSC
    if (_mesh->getProperties().existsPropertyVector<unsigned char>(
            "vtkGhostType", MeshLib::MeshItemType::Cell, 1))
    {
        auto* ghost_cell_property =
            _mesh->getProperties().getPropertyVector<unsigned char>(
                "vtkGhostType", MeshLib::MeshItemType::Cell, 1);
        if (ghost_cell_property)
        {
            const_cast<MeshLib::PropertyVector<unsigned char>*>(
                ghost_cell_property)
                ->is_for_output = true;
        }
    }
    else
    {
        DBUG("No vtkGhostType data in mesh '{}'.", _mesh->getName());
    }
#endif

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

#ifdef USE_PETSC
    if constexpr (std::is_same_v<UnstructuredGridWriter,
                                 vtkXMLPUnstructuredGridWriter>)
    {
        // Set the writer controller to same communicator as OGS
        vtkSmartPointer<vtkMPICommunicator> vtk_comm =
            vtkSmartPointer<vtkMPICommunicator>::New();
        MPI_Comm mpi_comm = MPI_COMM_WORLD;
        vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
        vtk_comm->InitializeExternal(&vtk_opaque_comm);

        vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl =
            vtkSmartPointer<vtkMPIController>::New();
        vtk_mpi_ctrl->SetCommunicator(vtk_comm);

        vtuWriter->SetController(vtk_mpi_ctrl);

        vtuWriter->SetGhostLevel(1);
        vtuWriter->SetNumberOfPieces(num_partitions);
        vtuWriter->SetStartPiece(rank);
        vtuWriter->SetEndPiece(rank);
    }
#endif

#ifdef VTK_USE_64BIT_IDS
    vtuWriter->SetHeaderTypeToUInt64();
    // set SetIdTypeToInt64() as well?
#endif

    return (vtuWriter->Write() > 0);
}

}  // end namespace IO
}  // end namespace MeshLib
