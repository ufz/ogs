/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-02
 * \brief  Implementation of the VtkPolylinesSource class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkPolylinesSource.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// GeoLib
#include "Color.h"
#include "Polyline.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

vtkStandardNewMacro(VtkPolylinesSource);
vtkCxxRevisionMacro(VtkPolylinesSource, "$Revision$");

VtkPolylinesSource::VtkPolylinesSource()
	: _polylines(NULL)
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	GetProperties()->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
}

VtkPolylinesSource::~VtkPolylinesSource()
{
}

void VtkPolylinesSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_polylines->size() == 0)
		return;

	for (std::vector<GeoLib::Polyline*>::const_iterator it = _polylines->begin();
	     it != _polylines->end(); ++it)
	{
		os << indent << "== Polyline ==" << "\n";
		int numPoints = (*it)->getNumberOfPoints();
		for (int i = 0; i < numPoints; i++)
		{
			const GeoLib::Point* point = (*it)->getPoint(i);
			const double* coords = point->getCoords();
			os << indent << "Point " << i << " (" << coords[0] << ", " << coords[1] <<
			", " << coords[2] << ")\n";
		}
	}
}

int VtkPolylinesSource::RequestData( vtkInformation* request,
                                     vtkInformationVector** inputVector,
                                     vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (!_polylines)
		return 0;
	if (_polylines->size() == 0)
	{
		ERR("VtkPolylineSource::RequestData(): Size of polyline vector is 0");
		return 0;
	}

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output =
	        vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();

	//newPoints->Allocate(numPoints);
	//newLines->Allocate(newLines->EstimateSize(numLines, 2));

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	vtkSmartPointer<vtkIntArray> plyIDs = vtkSmartPointer<vtkIntArray>::New();
	plyIDs->SetNumberOfComponents(1);
	plyIDs->SetName("PolylineIDs");

	int lastMaxIndex = 0;
	//for (std::vector<GeoLib::Polyline*>::const_iterator it = _polylines->begin();
	//	it != _polylines->end(); ++it)
	for (size_t j = 0; j < _polylines->size(); j++)
	{
		const int numPoints = (*_polylines)[j]->getNumberOfPoints();
		//const int numLines = numPoints - 1;

		// Generate points
		for (int i = 0; i < numPoints; i++)
		{
			const GeoLib::Point* point = (*_polylines)[j]->getPoint(i);
			const double* coords = point->getCoords();
			newPoints->InsertNextPoint(coords);
		}

		// Generate lines
		newLines->InsertNextCell(numPoints);
		plyIDs->InsertNextValue(j);
		for (int i = 0; i < numPoints; i++)
			newLines->InsertCellPoint(i + lastMaxIndex);

		lastMaxIndex += numPoints;
	}

	output->SetPoints(newPoints);
	output->SetLines(newLines);
	output->GetCellData()->AddArray(plyIDs);
	output->GetCellData()->SetActiveAttribute("PolylineIDs", vtkDataSetAttributes::SCALARS);
	output->Squeeze();

	return 1;
}

int VtkPolylinesSource::RequestInformation( vtkInformation* request,
                                            vtkInformationVector** inputVector,
                                            vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

	return 1;
}

void VtkPolylinesSource::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}
