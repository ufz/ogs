/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  VtkMappedMeshSource is a souce class to transform OGS meshes into complete
 *         vtkUnstructuredGrids.
 * Usage:
 * \code
 * vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
 * vtkSource->SetMesh(mesh);
 * vtkSource->Update();
 * vtkUnstructuredGrid* output = vtkSource->GetOutput();
 * \endcode
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    static VtkMappedMeshSource *New();
    vtkTypeMacro(VtkMappedMeshSource, vtkUnstructuredGridAlgorithm);
    void PrintSelf(std::ostream &os, vtkIndent indent) override;

    /// Sets the mesh. Calling is mandatory
    void SetMesh(const MeshLib::Mesh* mesh) { this->_mesh = mesh; this->Modified(); }

    /// Returns the mesh.
    const MeshLib::Mesh* GetMesh() const { return _mesh; }

protected:
    VtkMappedMeshSource();

    int ProcessRequest(vtkInformation *request, vtkInformationVector **inputVector,
                       vtkInformationVector *outputVector) override;
    int RequestData(vtkInformation *, vtkInformationVector **,
                    vtkInformationVector *) override;
    int RequestInformation(vtkInformation *, vtkInformationVector **,
                           vtkInformationVector *) override;

private:
    VtkMappedMeshSource(const VtkMappedMeshSource &); // Not implemented.
    void operator=(const VtkMappedMeshSource &);      // Not implemented.

    /// Adds a zero-copy vtk array wrapper.
    /// \param properties MeshLib::Properties object
    /// \param prop_name The name of the property vector to be mapped
    template <typename T>
    bool addProperty(MeshLib::Properties const& properties,
                     std::string const& prop_name) const
    {
        if (!properties.existsPropertyVector<T>(prop_name))
            return false;
        // TODO: Hack removing const
        auto* propertyVector = const_cast<MeshLib::PropertyVector<T> *>(
            properties.getPropertyVector<T>(prop_name));
        if(!propertyVector)
            return false;

        vtkNew<vtkAOSDataArrayTemplate<T> > dataArray;
        const bool hasArrayOwnership = false;
        dataArray->SetArray(propertyVector->data(),
            static_cast<vtkIdType>(propertyVector->size()),
            static_cast<int>(!hasArrayOwnership));
        dataArray->SetNumberOfComponents(propertyVector->getNumberOfComponents());
        dataArray->SetName(prop_name.c_str());

        if(propertyVector->getMeshItemType() == MeshLib::MeshItemType::Node)
            this->PointData->AddArray(dataArray.GetPointer());
        else if(propertyVector->getMeshItemType() == MeshLib::MeshItemType::Cell)
            this->CellData->AddArray(dataArray.GetPointer());
        else if(propertyVector->getMeshItemType() == MeshLib::MeshItemType::IntegrationPoint)
            this->FieldData->AddArray(dataArray.GetPointer());

        return true;
    }

    const MeshLib::Mesh* _mesh;

    int NumberOfDimensions;
    int NumberOfNodes;

    vtkNew<vtkPoints> Points;
    vtkNew<vtkPointData> PointData;
    vtkNew<vtkCellData> CellData;
    vtkNew<vtkFieldData> FieldData;
};

} // Namespace MeshLib
