/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-03
 * \brief  Implementation of the VtkSurfacesSource class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "Color.h"
#include "VtkSurfacesSource.h"
#include <limits>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkProperty.h>

vtkStandardNewMacro(VtkSurfacesSource);

VtkSurfacesSource::VtkSurfacesSource()
	: _surfaces(NULL)
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);
	//this->SetColorBySurface(true);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	vtkProperty* vtkProps = GetProperties();
	vtkProps->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	vtkProps->SetEdgeVisibility(0);
}

void VtkSurfacesSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_surfaces->size() == 0)
		return;

	os << indent << "== VtkSurfacesSource ==" << "\n";
}

int VtkSurfacesSource::RequestData( vtkInformation* request,
                                    vtkInformationVector** inputVector,
                                    vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	const int nSurfaces = _surfaces->size();
	if (nSurfaces == 0)
		return 0;

	const std::vector<GeoLib::Point*>* surfacePoints = (*_surfaces)[0]->getPointVec();
	size_t nPoints = surfacePoints->size();

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output =
	        vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	newPoints->SetNumberOfPoints(nPoints);

	vtkSmartPointer<vtkCellArray> newPolygons = vtkSmartPointer<vtkCellArray>::New();
	//newPolygons->Allocate(nSurfaces);

	vtkSmartPointer<vtkIntArray> sfcIDs = vtkSmartPointer<vtkIntArray>::New();
	sfcIDs->SetNumberOfComponents(1);
	sfcIDs->SetName("SurfaceIDs");

	for (size_t i = 0; i < nPoints; ++i)
	{
		const double* coords = const_cast<double*>((*surfacePoints)[i]->getCoords());
		newPoints->SetPoint(i, coords);
	}

	vtkIdType count(0);
	for (std::vector<GeoLib::Surface*>::const_iterator it = _surfaces->begin();
	     it != _surfaces->end(); ++it)
	{
		const size_t nTriangles = (*it)->getNTriangles();

		for (size_t i = 0; i < nTriangles; ++i)
		{
			vtkTriangle* new_tri = vtkTriangle::New();
			new_tri->GetPointIds()->SetNumberOfIds(3);

			const GeoLib::Triangle* triangle = (**it)[i];
			for (size_t j = 0; j < 3; ++j)
				new_tri->GetPointIds()->SetId(j, ((*triangle)[j]));
			newPolygons->InsertNextCell(new_tri);
			sfcIDs->InsertNextValue(count);
			new_tri->Delete();
		}
		count++;
	}

	output->SetPoints(newPoints);
	output->SetPolys(newPolygons);
	output->GetCellData()->AddArray(sfcIDs);
	output->GetCellData()->SetActiveAttribute("SurfaceIDs", vtkDataSetAttributes::SCALARS);
	output->Squeeze();

	return 1;
}

int VtkSurfacesSource::RequestInformation( vtkInformation* /*request*/,
                                           vtkInformationVector** /*inputVector*/,
                                           vtkInformationVector* /*outputVector*/ )
{
	return 1;
}

void VtkSurfacesSource::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);
	(*_algorithmUserProperties)[name] = value;
}

