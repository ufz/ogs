/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  Implementation of the VtkMappedMeshSource class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "VtkMappedMeshSource.h"

#include <vector>

#include <vtkCellType.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

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
    : NumberOfDimensions(0), NumberOfNodes(0)
{
    this->SetNumberOfInputPorts(0);
}

int VtkMappedMeshSource::ProcessRequest(vtkInformation* request,
                                        vtkInformationVector** inputVector,
                                        vtkInformationVector* outputVector)
{
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
        return this->RequestData(request, inputVector, outputVector);

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
        return this->RequestInformation(request, inputVector, outputVector);

    return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

int VtkMappedMeshSource::RequestData(vtkInformation*,
                                     vtkInformationVector**,
                                     vtkInformationVector* outputVector)
{
    vtkSmartPointer<vtkInformation> outInfo =
        outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkUnstructuredGrid> output =
        vtkUnstructuredGrid::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));

    if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) >
        0)
        return 1;

    // Points
    this->Points->Reset();

    vtkNew<VtkMeshNodalCoordinatesTemplate<double>> nodeCoords;
    nodeCoords->SetNodes(_mesh->getNodes());
    this->Points->SetData(nodeCoords.GetPointer());
    output->SetPoints(this->Points.GetPointer());

    // Cells
    auto elems = _mesh->getElements();
    output->Allocate(elems.size());
    for (std::size_t cell = 0; cell < elems.size(); cell++)
    {
        auto cellType = OGSToVtkCellType(elems[cell]->getCellType());

        const MeshLib::Element* const elem = (elems)[cell];
        const unsigned numNodes(elem->getNumberOfNodes());
        const MeshLib::Node* const* nodes = elems[cell]->getNodes();
        vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
        ptIds->SetNumberOfIds(numNodes);

        for (unsigned i = 0; i < numNodes; ++i)
            ptIds->SetId(i, nodes[i]->getID());

        if (cellType == VTK_WEDGE)
        {
            for (unsigned i = 0; i < 3; ++i)
            {
                const unsigned prism_swap_id = ptIds->GetId(i);
                ptIds->SetId(i, ptIds->GetId(i + 3));
                ptIds->SetId(i + 3, prism_swap_id);
            }
        }
        else if (cellType == VTK_QUADRATIC_WEDGE)
        {
            std::array<vtkIdType, 15> ogs_nodeIds;
            for (unsigned i = 0; i < 15; ++i)
                ogs_nodeIds[i] = ptIds->GetId(i);
            for (unsigned i = 0; i < 3; ++i)
            {
                ptIds->SetId(i, ogs_nodeIds[i + 3]);
                ptIds->SetId(i + 3, ogs_nodeIds[i]);
            }
            for (unsigned i = 0; i < 3; ++i)
                ptIds->SetId(6 + i, ogs_nodeIds[8 - i]);
            for (unsigned i = 0; i < 3; ++i)
                ptIds->SetId(9 + i, ogs_nodeIds[14 - i]);
            ptIds->SetId(12, ogs_nodeIds[9]);
            ptIds->SetId(13, ogs_nodeIds[11]);
            ptIds->SetId(14, ogs_nodeIds[10]);
        }

        output->InsertNextCell(cellType, ptIds);
    }

    // Arrays
    MeshLib::Properties const& properties = _mesh->getProperties();
    std::vector<std::string> const& propertyNames =
        properties.getPropertyVectorNames();

    for (std::vector<std::string>::const_iterator name = propertyNames.cbegin();
         name != propertyNames.cend(); ++name)
    {
        if (addProperty<double>(properties, *name))
            continue;
        if (addProperty<int>(properties, *name))
            continue;
        if (addProperty<unsigned>(properties, *name))
            continue;
        if (addProperty<std::size_t>(properties, *name))
            continue;
        if (addProperty<char>(properties, *name))
            continue;

        DBUG("Mesh property \"%s\" with unknown data type.", *name->c_str());
    }

    output->GetPointData()->ShallowCopy(this->PointData.GetPointer());
    output->GetCellData()->ShallowCopy(this->CellData.GetPointer());
    output->GetFieldData()->ShallowCopy(this->FieldData.GetPointer());
    return 1;
}

int VtkMappedMeshSource::RequestInformation(vtkInformation*,
                                            vtkInformationVector**,
                                            vtkInformationVector*)
{
    this->NumberOfDimensions = 3;
    this->NumberOfNodes = _mesh->getNumberOfNodes();

    return 1;
}

}  // Namespace MeshLib
