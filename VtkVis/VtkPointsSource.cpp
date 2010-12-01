/**
 * \file VtkPointsSource.cpp
 * 3/2/2010 LB Initial implementation
 *
 * Implementation of VtkPointsSource
 */

// ** INCLUDES **
#include "VtkPointsSource.h"

#include <vtkCellArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

vtkStandardNewMacro(VtkPointsSource);
vtkCxxRevisionMacro(VtkPointsSource, "$Revision$");

VtkPointsSource::VtkPointsSource()
: _points(NULL)
{
	this->SetNumberOfInputPorts(0);

	const GEOLIB::Color* c = GEOLIB::getRandomColor();
	GetProperties()->SetColor((*c)[0]/255.0,(*c)[1]/255.0,(*c)[2]/255.0);
}

void VtkPointsSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_points->size() == 0)
		return;

	os << indent << "== VtkPointsSource ==" << "\n";

	int i = 0;
	for (std::vector<GEOLIB::Point*>::const_iterator it = _points->begin();
		it != _points->end(); ++it)
	{
		const double* coords = (*it)->getData();
		os << indent << "Point " << i <<" (" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")\n";
		i++;
	}
}

int VtkPointsSource::RequestData( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (!_points)
		return 0;
	int numPoints = _points->size();
	if (numPoints == 0)
	{
		std::cout << "ERROR in VtkPointsSource::RequestData : Size of point vector is 0" << std::endl;
		return 0;
	}

	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkPoints* newPoints = vtkPoints::New();
	vtkCellArray* newVerts = vtkCellArray::New();
	newPoints->Allocate(numPoints);
	newVerts->Allocate(numPoints);

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	// Generate points and vertices
	for (std::vector<GEOLIB::Point*>::const_iterator it = _points->begin();
		it != _points->end(); ++it)
	{
		double coords[3] = {(*(*it))[0], (*(*it))[1], (*(*it))[2]};
		vtkIdType pid[1];
		pid[0] = newPoints->InsertNextPoint(coords);
		newVerts->InsertNextCell(1, pid);
	}

	output->SetPoints(newPoints);
	newPoints->Delete();

	output->SetVerts(newVerts);
	newVerts->Delete();

	return 1;
}

int VtkPointsSource::RequestInformation( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

	return 1;
}

void VtkPointsSource::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}
