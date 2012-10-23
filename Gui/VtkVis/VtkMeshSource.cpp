/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VtkMeshSource.cpp
 *
 * Created on 2010-03-19 by Karsten Rink
 */

#include "VtkColorLookupTable.h"
#include "VtkMeshSource.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

#include "Color.h"

// ** VTK INCLUDES **
#include "vtkObjectFactory.h"
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProperty.h>

// OGS Cell Types
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkWedge.h> // == Prism
#include <vtkPyramid.h>

vtkStandardNewMacro(VtkMeshSource);
vtkCxxRevisionMacro(VtkMeshSource, "$Revision$");

VtkMeshSource::VtkMeshSource() : _matName("MaterialIDs")
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	vtkProperty* vtkProps = GetProperties();
	vtkProps->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	vtkProps->SetEdgeVisibility(1);
}

VtkMeshSource::~VtkMeshSource()
{
}

void VtkMeshSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_grid == NULL)
		return;
	const std::vector<MeshLib::Node*> nodes = _grid->getNodes();
	const std::vector<MeshLib::Element*> elems = _grid->getElements();
	if (nodes.empty() || elems.empty() )
		return;

	os << indent << "== VtkMeshSource ==" << "\n";

	int i = 0;
	for (std::vector<MeshLib::Node*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		os << indent << "Point " << i << " (" << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << ")" << std::endl;
	}

	i = 0;
	for (std::vector<MeshLib::Element*>::const_iterator it = elems.begin(); it != elems.end(); ++it)
	{
		os << indent << "Element " << i << ": ";
		for (unsigned t = 0; t < (*it)->getNNodes(); ++t)
			os << (*it)->getNode(t)->getID() << " ";
		os << std::endl;
	}
}

int VtkMeshSource::RequestData( vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (_grid == NULL)
		return 0;
	const std::vector<MeshLib::Node*> nodes = _grid->getNodes();
	const std::vector<MeshLib::Element*> elems = _grid->getElements();

	const size_t nPoints = _grid->getNNodes();
	const size_t nElems  = _grid->getNElements();
	if (nPoints == 0 || nElems == 0)
		return 0;

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkUnstructuredGrid> output = vtkUnstructuredGrid::SafeDownCast(
	        outInfo->Get(vtkDataObject::DATA_OBJECT()));
	output->Allocate(nElems);

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	// Insert grid points
	vtkSmartPointer<vtkPoints> gridPoints = vtkSmartPointer<vtkPoints>::New();
	gridPoints->Allocate(nPoints);
	// Generate mesh nodes
	for (unsigned i = 0; i < nPoints; ++i)
		gridPoints->InsertPoint(i, (*nodes[i])[0], (*nodes[i])[1], (*nodes[i])[2]);

	// Generate attribute vector for material groups
	vtkSmartPointer<vtkIntArray> materialIDs = vtkSmartPointer<vtkIntArray>::New();
	materialIDs->SetName(_matName);
	materialIDs->SetNumberOfComponents(1);
	materialIDs->SetNumberOfTuples(nElems);

	// Generate mesh elements
	for (unsigned i = 0; i < nElems; ++i)
	{
		int type(0);
		const MeshLib::Element* elem (elems[i]);

		switch (elem->getGeomType())
		{
		case MshElemType::EDGE:
			type = 3;
			break;
		case MshElemType::TRIANGLE:
			type = 5;
			break;
		case MshElemType::QUAD:
			type = 9;
			break;
		case MshElemType::HEXAHEDRON:
			type = 12;
			break;
		case MshElemType::TETRAHEDRON:
			type = 10;
			break;
		case MshElemType::PRISM:
			type = 13;
			break;
		case MshElemType::PYRAMID:
			type = 14;
			break;
		default: // if none of the above can be applied
			std::cout << "Error in VtkMeshSource::RequestData() - Unknown element type " << MshElemType2String(elem->getGeomType()) << "." << std::endl;
			return 0;
		}

		materialIDs->InsertValue(i, elem->getValue());
		vtkIdList* point_ids = vtkIdList::New();

		const unsigned nElemNodes (elem->getNNodes());
		for (unsigned j = 0; j < nElemNodes; ++j)
			point_ids->InsertNextId(elem->getNode(j)->getID());

		output->InsertNextCell(type, point_ids);
	}

	output->SetPoints(gridPoints);

	output->GetCellData()->AddArray(materialIDs);
	output->GetCellData()->SetActiveAttribute(_matName, vtkDataSetAttributes::SCALARS);

	return 1;
}

void VtkMeshSource::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	(*_algorithmUserProperties)[name] = value;
}
