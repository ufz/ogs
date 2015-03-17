/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  Implementation of the VtkMappedMeshSource class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMappedMeshSource.h"

#include "logog/include/logog.hpp"

#include <vtkCellData.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>

#include "VtkMappedMesh.h"
#include "VtkMeshNodalCoordinatesTemplate.h"
#include "VtkMappedPropertyVectorTemplate.h"

namespace InSituLib {

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

VtkMappedMeshSource::~VtkMappedMeshSource()
{

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

	// Elements
	vtkNew<VtkMappedMesh> elems;
	elems->GetImplementation()->SetNodes(_mesh->getNodes());
	elems->GetImplementation()->SetElements(_mesh->getElements());
	elems->SetPoints(this->Points.GetPointer());
	output->Allocate(elems->GetNumberOfCells());
	output->ShallowCopy(elems.GetPointer());

	// Arrays
	MeshLib::Properties const & properties = _mesh->getProperties();
	std::vector<std::string> const& propertyNames = properties.getPropertyVectorNames();

	for(std::vector<std::string>::const_iterator it = propertyNames.cbegin();
		it != propertyNames.cend(); ++it)
	{
		if (addProperty<double>(properties, *it))
			continue;

		if (addProperty<int>(properties, *it))
			continue;

		DBUG ("Mesh property \"%s\" with unknown data type.", *it->c_str());
	}

	output->GetPointData()->ShallowCopy(this->PointData.GetPointer());
	output->GetCellData()->ShallowCopy(this->CellData.GetPointer());

	return 1;
}

template<typename T> bool VtkMappedMeshSource::addProperty(MeshLib::Properties const& properties,
                                                           std::string const& prop_name) const
{
	boost::optional<MeshLib::PropertyVector<T> const &>
		propertyVector(properties.getPropertyVector<T>(prop_name));
	if(!propertyVector)
		return false;

	vtkNew<VtkMappedPropertyVectorTemplate<T> > dataArray;
	dataArray->SetPropertyVector(const_cast<MeshLib::PropertyVector<T> &>(*propertyVector));
	dataArray->SetName(prop_name.c_str());

	if(propertyVector->getMeshItemType() == MeshLib::MeshItemType::Node)
		this->PointData->AddArray(dataArray.GetPointer());
	else if(propertyVector->getMeshItemType() == MeshLib::MeshItemType::Cell)
		this->CellData->AddArray(dataArray.GetPointer());

	return true;
}

int VtkMappedMeshSource::RequestInformation(
	vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
	this->NumberOfDimensions = 3;
	this->NumberOfNodes = _mesh->getNNodes();

	return 1;
}

} // Namespace InSituLib
