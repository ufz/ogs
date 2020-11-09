/**
 * \file
 * \author Tobias Meisel
 * \date   2020-10-06
 * \brief  Implementation of WriteXdmf3 function.
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
bool writeXdmf3(const MeshLib::Mesh& mesh,
                std::filesystem::path const& file_path)
{
    vtkSmartPointer<vtkXdmf3Writer> writer =
        vtkSmartPointer<vtkXdmf3Writer>::New();

    writer->SetFileName(file_path.string().c_str());
    vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
    vtkSource->SetMesh(&mesh);
    vtkSource->Update();
    writer->SetInputData(vtkSource->GetOutput());
    writer->Write();
    return true;
}
}  // end namespace IO
}  // end namespace MeshLib