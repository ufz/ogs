/**
 * \file
 * \author Tobias Meisel
 * \date   2020-10-06
 * \brief  Implementation of WriteXdmf function.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeXdmf.h"

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkXdmf3Writer.h>

#include "MeshLib/Vtk/VtkMappedMeshSource.h"

namespace MeshLib
{
namespace IO
{
bool writeXdmf3(const MeshLib::Mesh& mesh, std::string const &file_name)
{
        #ifdef USE_PETSC
            int rank;
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            int mpi_size;
            MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
            vtkSmartPointer<vtkPXdmf3Writer> writer = vtkSmartPointer<vtkPXdmf3Writer>::New(); // open file handle

        #else
            vtkSmartPointer<vtkXdmf3Writer> writer = vtkSmartPointer<vtkXdmf3Writer>::New(); // open file handle
        #endif

            writer->SetFileName(file_name.c_str());
            vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
            vtkSource->SetMesh(&mesh);
            vtkSource->Update();
            writer->SetInputData(vtkSource->GetOutput());

        #ifdef USE_PETSC
            writer->SetGhostLevel(1);
            writer->SetNumberOfPieces(num_partitions);
            writer->SetStartPiece(rank);
            writer->SetEndPiece(rank);
        #endif
            writer->Write();
                // close file handle
        return 0;
}
} //end namespace IO
} //end namespace MeshLib