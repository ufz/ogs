/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  Implementation of the VtkMappedMeshSource class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "VtkMappedMeshSource.h"

#include <vtkCellType.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/VtkOGSEnum.h"
#include "VtkMeshNodalCoordinatesTemplate.h"

namespace MeshLib
{
vtkStandardNewMacro(VtkMappedMeshSource)

    void VtkMappedMeshSource::PrintSelf(std::ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Mesh: " << (_mesh ? _mesh->getName() : "(none)") << endl;
}

VtkMappedMeshSource::VtkMappedMeshSource()

{
    this->SetNumberOfInputPorts(0);
}

int VtkMappedMeshSource::ProcessRequest(vtkInformation* request,
                                        vtkInformationVector** inputVector,
                                        vtkInformationVector* outputVector)
{
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
        return this->RequestData(request, inputVector, outputVector);
    }

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
        return this->RequestInformation(request, inputVector, outputVector);
    }

    return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

int VtkMappedMeshSource::RequestData(vtkInformation* /*request*/,
                                     vtkInformationVector** /*inputVector*/,
                                     vtkInformationVector* outputVector)
{
    vtkSmartPointer<vtkInformation> outInfo =
        outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkUnstructuredGrid> output =
        vtkUnstructuredGrid::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));

    if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) >
        0)
    {
        return 1;
    }

    // Points
    this->Points->Reset();

    vtkSmartPointer<VtkMeshNodalCoordinatesTemplate<double>> nodeCoords =
        vtkSmartPointer<VtkMeshNodalCoordinatesTemplate<double>>::New();
    nodeCoords->SetNodes(_mesh->getNodes());
    this->Points->SetData(nodeCoords);
    output->SetPoints(this->Points.GetPointer());

    // Cells
    auto elems = _mesh->getElements();
    output->Allocate(elems.size());
    for (auto& cell : elems)
    {
        auto cellType = OGSToVtkCellType(cell->getCellType());

        const MeshLib::Element* const elem = cell;
        const unsigned numNodes(elem->getNumberOfNodes());
        const MeshLib::Node* const* nodes = cell->getNodes();
        vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
        ptIds->SetNumberOfIds(numNodes);

        for (unsigned i = 0; i < numNodes; ++i)
        {
            ptIds->SetId(i, nodes[i]->getID());
        }

        if (cellType == VTK_WEDGE)
        {
            for (unsigned i = 0; i < 3; ++i)
            {
                const auto prism_swap_id = ptIds->GetId(i);
                ptIds->SetId(i, ptIds->GetId(i + 3));
                ptIds->SetId(i + 3, prism_swap_id);
            }
        }

        output->InsertNextCell(cellType, ptIds);
    }

    // Arrays
    for (auto [name, property] : _mesh->getProperties())
    {
        if (auto p = dynamic_cast<PropertyVector<double>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<float>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<int>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<unsigned>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<long>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<long long>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p =
                     dynamic_cast<PropertyVector<unsigned long>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<unsigned long long>*>(
                     property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<std::size_t>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<char>*>(property))
        {
            addProperty(*p);
        }
        else if (auto p =
                     dynamic_cast<PropertyVector<unsigned char>*>(property))
        {
            addProperty(*p);
        }
        else
        {
            OGS_FATAL(
                "Mesh property '{:s}' of unhandled data type '{:s}'. Please "
                "check the data type of the mesh properties. The available "
                "data types are:"
                "\n\t double,"
                "\n\t float,"
                "\n\t int,"
                "\n\t unsigned,"
                "\n\t long,"
                "\n\t long long,"
                "\n\t unsigned long long,"
                "\n\t char,",
                "\n\t unsigned char.",
                property->getPropertyName(),
                typeid(*property).name());
        }
    }

    output->GetPointData()->ShallowCopy(this->PointData.GetPointer());
    output->GetCellData()->ShallowCopy(this->CellData.GetPointer());
    output->GetFieldData()->ShallowCopy(this->FieldData.GetPointer());
    return 1;
}

int VtkMappedMeshSource::RequestInformation(
    vtkInformation* /*request*/,
    vtkInformationVector** /*inputVector*/,
    vtkInformationVector* /*outputVector*/)
{
    this->NumberOfDimensions = 3;
    this->NumberOfNodes = _mesh->getNumberOfNodes();

    return 1;
}

template <typename T>
void VtkMappedMeshSource::addProperty(
    MeshLib::PropertyVector<T> const& property) const
{
    vtkNew<vtkAOSDataArrayTemplate<T>> dataArray;
    const bool hasArrayOwnership = false;
    dataArray->SetArray(const_cast<T*>(property.data()),
                        static_cast<vtkIdType>(property.size()),
                        static_cast<int>(!hasArrayOwnership));
    dataArray->SetNumberOfComponents(property.getNumberOfGlobalComponents());
    dataArray->SetName(property.getPropertyName().c_str());

    if (property.getMeshItemType() == MeshLib::MeshItemType::Node)
    {
        this->PointData->AddArray(dataArray.GetPointer());
    }
    else if (property.getMeshItemType() == MeshLib::MeshItemType::Cell)
    {
        this->CellData->AddArray(dataArray.GetPointer());
    }
    else if (property.getMeshItemType() ==
             MeshLib::MeshItemType::IntegrationPoint)
    {
        this->FieldData->AddArray(dataArray.GetPointer());
    }
}
}  // Namespace MeshLib
