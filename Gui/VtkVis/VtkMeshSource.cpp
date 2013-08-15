/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-19
 * \brief  Implementation of the VtkMeshSource class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMeshSource.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Elements/Element.h"
#include "Mesh.h"
#include "Node.h"
#include "VtkColorLookupTable.h"

#include "Color.h"

// ** VTK INCLUDES **
#include "vtkObjectFactory.h"
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

// OGS Cell Types
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkPyramid.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkWedge.h> // == Prism

vtkStandardNewMacro(VtkMeshSource);
vtkCxxRevisionMacro(VtkMeshSource, "$Revision$");

VtkMeshSource::VtkMeshSource() :
		_grid(nullptr), _matName("MaterialIDs")
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	vtkProperty* vtkProps = GetProperties();
	vtkProps->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	delete c;
	vtkProps->SetEdgeVisibility(1);
}

VtkMeshSource::~VtkMeshSource()
{
}

void VtkMeshSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_grid == nullptr)
		return;
	const std::vector<MeshLib::Node*> nodes = _grid->getNodes();
	const std::vector<MeshLib::Element*> elems = _grid->getElements();
	if (nodes.empty() || elems.empty() )
		return;

	os << indent << "== VtkMeshSource ==" << "\n";

	int i = 0;
	for (std::vector<MeshLib::Node*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		os << indent << "Point " << i << " (" << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << ")\n";
	}

	i = 0;
	for (std::vector<MeshLib::Element*>::const_iterator it = elems.begin(); it != elems.end(); ++it)
	{
		os << indent << "Element " << i << ": ";
		for (unsigned t = 0; t < (*it)->getNNodes(); ++t)
			os << (*it)->getNode(t)->getID() << " ";
		os << "\n";
	}
}

int VtkMeshSource::RequestData( vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (_grid == nullptr)
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
		case MeshElemType::EDGE:
			type = 3;
			break;
		case MeshElemType::TRIANGLE:
			type = 5;
			break;
		case MeshElemType::QUAD:
			type = 9;
			break;
		case MeshElemType::HEXAHEDRON:
			type = 12;
			break;
		case MeshElemType::TETRAHEDRON:
			type = 10;
			break;
		case MeshElemType::PRISM:
			type = 13;
			break;
		case MeshElemType::PYRAMID:
			type = 14;
			break;
		default: // if none of the above can be applied
			ERR("VtkMeshSource::RequestData(): Unknown element type \"%s\".",
					MeshElemType2String(elem->getGeomType()).c_str());
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
