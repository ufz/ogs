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

#include <vtkDemandDrivenPipeline.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include "VtkMappedMesh.h"
#include "VtkMeshNodalCoordinatesTemplate.h"

namespace MeshLib {

vtkStandardNewMacro(VtkMappedMeshSource)

void VtkMappedMeshSource::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Mesh: " << (_mesh ? _mesh->getName() : "(none)") << endl;
}

VtkMappedMeshSource::VtkMappedMeshSource()
 : NumberOfDimensions(0),
   NumberOfNodes(0)
{
    this->SetNumberOfInputPorts(0);
}

int VtkMappedMeshSource::ProcessRequest(
    vtkInformation *request, vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
        return this->RequestData(request, inputVector, outputVector);

    if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
        return this->RequestInformation(request, inputVector, outputVector);

    return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

int VtkMappedMeshSource::RequestData(vtkInformation *,
                                     vtkInformationVector **,
                                     vtkInformationVector *outputVector)
{
    vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkUnstructuredGrid> output = vtkUnstructuredGrid::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
        return 1;

    // Points
    this->Points->Reset();

    vtkNew<VtkMeshNodalCoordinatesTemplate<double> > nodeCoords;
    nodeCoords->SetNodes(_mesh->getNodes());
    this->Points->SetData(nodeCoords.GetPointer());
    // output->SetPoints(this->Points.GetPointer()); // TODO: not necessary?

    // Elements
    vtkNew<VtkMappedMesh> elems;
    elems->GetImplementation()->SetNodes(_mesh->getNodes());
    elems->GetImplementation()->SetElements(_mesh->getElements());

     // Use the mapped point container for the block points
    elems->SetPoints(this->Points.GetPointer());

    output->Allocate(elems->GetNumberOfCells());
    output->ShallowCopy(elems.GetPointer());

    // Arrays
    MeshLib::Properties const & properties = _mesh->getProperties();
    std::vector<std::string> const& propertyNames = properties.getPropertyVectorNames();

    for(std::vector<std::string>::const_iterator name = propertyNames.cbegin(); name != propertyNames.cend(); ++name)
    {
        if (addProperty<double>(properties, *name)) continue;
        if (addProperty<int>(properties, *name)) continue;
        if (addProperty<unsigned>(properties, *name)) continue;
        if (addProperty<std::size_t>(properties, *name)) continue;
        if (addProperty<char>(properties, *name)) continue;

        DBUG ("Mesh property \"%s\" with unknown data type.", *name->c_str());
    }

    output->GetPointData()->ShallowCopy(this->PointData.GetPointer());
    output->GetCellData()->ShallowCopy(this->CellData.GetPointer());
    return 1;
}

int VtkMappedMeshSource::RequestInformation(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
    this->NumberOfDimensions = 3;
    this->NumberOfNodes = _mesh->getNumberOfNodes();

    return 1;
}

} // Namespace MeshLib
