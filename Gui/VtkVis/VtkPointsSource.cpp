/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-03
 * \brief  Implementation of the VtkPointsSource class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
// GeoLib
#include "Color.h"

#include "VtkPointsSource.h"

#include <vtkCellArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellData.h>
#include <vtkProperty.h>

vtkStandardNewMacro(VtkPointsSource);
vtkCxxRevisionMacro(VtkPointsSource, "$Revision$");

VtkPointsSource::VtkPointsSource()
	: _points(NULL)
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	GetProperties()->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	delete c;
}

void VtkPointsSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_points->size() == 0)
		return;

	os << indent << "== VtkPointsSource ==" << "\n";

	int i = 0;
	for (std::vector<GeoLib::Point*>::const_iterator it = _points->begin();
	     it != _points->end(); ++it)
	{
		const double* coords = (*it)->getCoords();
		os << indent << "Point " << i << " (" << coords[0] << ", " << coords[1] << ", " <<
		coords[2] << ")\n";
		i++;
	}
}

int VtkPointsSource::RequestData( vtkInformation* request,
                                  vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (!_points)
		return 0;
	int numPoints = _points->size();
	if (numPoints == 0)
	{
		std::cout << "ERROR in VtkPointsSource::RequestData : Size of point vector is 0" <<
		std::endl;
		return 0;
	}

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output =
	        vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
	newPoints->Allocate(numPoints);
	newVerts->Allocate(numPoints);

	vtkSmartPointer<vtkIntArray> pointIDs = vtkSmartPointer<vtkIntArray>::New();
	pointIDs->SetNumberOfComponents(1);
	pointIDs->SetName("PointIDs");

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	// Generate points and vertices
	int i = 0;
	for (std::vector<GeoLib::Point*>::const_iterator it = _points->begin();
	     it != _points->end(); ++it)
	{
		double coords[3] = {(*(*it))[0], (*(*it))[1], (*(*it))[2]};
		vtkIdType pid = newPoints->InsertNextPoint(coords);
		newVerts->InsertNextCell(1, &pid);

		pointIDs->InsertNextValue(i);
		i++;
	}

	output->SetPoints(newPoints);
	output->SetVerts(newVerts);
	output->GetCellData()->AddArray(pointIDs);
	output->GetCellData()->SetActiveAttribute("PointIDs", vtkDataSetAttributes::SCALARS);

	return 1;
}

int VtkPointsSource::RequestInformation( vtkInformation* request,
                                         vtkInformationVector** inputVector,
                                         vtkInformationVector* outputVector )
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
