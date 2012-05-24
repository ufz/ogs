/**
 * \file VtkMeshSource.cpp
 * 19/03/2010 KR Initial implementation
 *
 */

#include "VtkColorLookupTable.h"
#include "VtkMeshSource.h"

// ** VTK INCLUDES **
#include "vtkObjectFactory.h"
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnsignedCharArray.h>
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

	const GEOLIB::Color* c = GEOLIB::getRandomColor();
	vtkProperty* vtkProps = GetProperties();
	vtkProps->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	vtkProps->SetEdgeVisibility(1);
}

VtkMeshSource::~VtkMeshSource()
{
	delete _grid;
}

void VtkMeshSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_grid == NULL)
		return;
	const std::vector<GEOLIB::Point*>* nodes = _grid->getNodes();
	const std::vector<GridAdapter::Element*>* elems = _grid->getElements();
	if (nodes->empty() || elems->empty() )
		return;

	os << indent << "== VtkMeshSource ==" << "\n";

	int i = 0;
	for (std::vector<GEOLIB::Point*>::const_iterator it = nodes->begin();
	     it != nodes->end(); ++it)
	{
		os << indent << "Point " << i << " (" << (*it)[0] << ", " << (*it)[1] << ", " <<
		(*it)[2] << ")" << std::endl;
		i++;
	}

	i = 0;
	for (std::vector<GridAdapter::Element*>::const_iterator it = elems->begin();
	     it != elems->end(); ++it)
	{
		os << indent << "Element " << i << ": ";
		for (size_t t = 0; t < (*it)->nodes.size(); t++)
			os << (*it)->nodes[t] << " ";
		os << std::endl;
		i++;
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
	const std::vector<GEOLIB::Point*>* nodes = _grid->getNodes();
	const std::vector<GridAdapter::Element*>* elems = _grid->getElements();

	const size_t nPoints = nodes->size();
	const size_t nElems  = elems->size();
	//size_t nElemNodes = 0;
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
	for (size_t i = 0; i < nPoints; i++)
		gridPoints->InsertPoint(i, (*(*nodes)[i])[0], (*(*nodes)[i])[1], (*(*nodes)[i])[2]);

	// Generate attribute vector for material groups
	vtkSmartPointer<vtkIntArray> materialIDs = vtkSmartPointer<vtkIntArray>::New();
	materialIDs->SetName(_matName);
	materialIDs->SetNumberOfComponents(1);
	materialIDs->SetNumberOfTuples(nElems);

	// Generate mesh elements
	for (size_t i = 0; i < nElems; i++)
	{
		vtkCell* newCell;
		const GridAdapter::Element* elem = (*elems)[i];

		switch (elem->type)
		{
		case MshElemType::LINE:
			newCell = vtkLine::New();
			break;
		case MshElemType::TRIANGLE:
			newCell = vtkTriangle::New();
			break;
		case MshElemType::QUAD:
			newCell = vtkQuad::New();
			break;
		case MshElemType::HEXAHEDRON:
			newCell = vtkHexahedron::New();
			break;
		case MshElemType::TETRAHEDRON:
			newCell = vtkTetra::New();
			break;
		case MshElemType::PRISM:
			newCell = vtkWedge::New();
			break;
		case MshElemType::PYRAMID:
			newCell = vtkPyramid::New();
			break;
		default: // if none of the above can be applied
			std::cout << "Error in VtkMeshSource::RequestData() - Unknown element type " << MshElemType2String(elem->type) << "." << std::endl;
			return 0;
		}

		materialIDs->InsertValue(i,(elem->material));

		const size_t nElemNodes (elem->nodes.size());
		for (size_t j = 0; j < nElemNodes; j++)
			newCell->GetPointIds()->SetId(j, elem->nodes[nElemNodes - 1 - j]);

		output->InsertNextCell(newCell->GetCellType(), newCell->GetPointIds());
		newCell->Delete();
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
