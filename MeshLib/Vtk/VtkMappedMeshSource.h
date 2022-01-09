/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  VtkMappedMeshSource is a source class to transform OGS meshes into
 * complete vtkUnstructuredGrids. Usage: \code
 * vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
 * vtkSource->SetMesh(mesh);
 * vtkSource->Update();
 * vtkUnstructuredGrid* output = vtkSource->GetOutput();
 * \endcode
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>

#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"



namespace MeshLib {

class Mesh;

/// Adapter which maps a MeshLib::Mesh to a vtkUnstructuredGridAlgorithm.
/// Allows for zero-copy access of the mesh from the visualization side.
class VtkMappedMeshSource final : public vtkUnstructuredGridAlgorithm
{
public:
    static VtkMappedMeshSource* New();
    vtkTypeMacro(VtkMappedMeshSource, vtkUnstructuredGridAlgorithm);
    void PrintSelf(std::ostream& os, vtkIndent indent) override;

    VtkMappedMeshSource(const VtkMappedMeshSource&) = delete;
    void operator=(const VtkMappedMeshSource&) = delete;

    /// Sets the mesh. Calling is mandatory
    void SetMesh(const MeshLib::Mesh* mesh)
    {
        this->_mesh = mesh;
        this->Modified();
    }

    /// Returns the mesh.
    const MeshLib::Mesh* GetMesh() const { return _mesh; }

protected:
    VtkMappedMeshSource();

    int ProcessRequest(vtkInformation* request,
                       vtkInformationVector** inputVector,
                       vtkInformationVector* outputVector) override;
    int RequestData(vtkInformation* /*request*/,
                    vtkInformationVector** /*inputVector*/,
                    vtkInformationVector* /*outputVector*/) override;
    int RequestInformation(vtkInformation* /*request*/,
                           vtkInformationVector** /*inputVector*/,
                           vtkInformationVector* /*outputVector*/) override;

private:
    /// Adds a zero-copy vtk array wrapper.
    template <typename T>
    void addProperty(MeshLib::PropertyVector<T> const& property) const;

    const MeshLib::Mesh* _mesh{};

    int NumberOfDimensions{0};
    int NumberOfNodes{0};

    vtkNew<vtkPoints> Points;
    vtkNew<vtkPointData> PointData;
    vtkNew<vtkCellData> CellData;
    vtkNew<vtkFieldData> FieldData;
};

} // Namespace MeshLib
