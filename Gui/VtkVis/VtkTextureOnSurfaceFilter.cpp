/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-03
 * \brief  Implementation of the VtkTextureOnSurfaceFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkTextureOnSurfaceFilter.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkImageShiftScale.h>
#include <vtkImageAlgorithm.h>
#include <vtkProperty.h>
#include <vtkTexture.h>

#include "MathTools.h"
#include "VtkVisHelper.h"

vtkStandardNewMacro(VtkTextureOnSurfaceFilter);
vtkCxxRevisionMacro(VtkTextureOnSurfaceFilter, "$Revision$");

VtkTextureOnSurfaceFilter::VtkTextureOnSurfaceFilter() : _origin(0,0), _scalingFactor(0)
{
}

VtkTextureOnSurfaceFilter::~VtkTextureOnSurfaceFilter()
{
}

void VtkTextureOnSurfaceFilter::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);
}

int VtkTextureOnSurfaceFilter::RequestData( vtkInformation* request,
                                            vtkInformationVector** inputVector,
                                            vtkInformationVector* outputVector )
{
	(void)request;

	if (this->GetTexture() == nullptr)
	{
		ERR("VtkTextureOnSurfaceFilter::RequestData() - No texture specified.");
		return 0;
	}

	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	int dims[3];
	this->GetTexture()->GetInput()->GetDimensions(dims);
	size_t imgWidth = dims[0];
	size_t imgHeight = dims[1];

	std::pair<int, int> min((int)_origin.first, (int)_origin.second);
	std::pair<int, int> max((int)(_origin.first + (imgWidth * _scalingFactor)),
	                        (int)(_origin.second + (imgHeight * _scalingFactor)));

	//calculate texture coordinates
	vtkPoints* points = input->GetPoints();
	vtkSmartPointer<vtkFloatArray> textureCoordinates = vtkSmartPointer<vtkFloatArray>::New();
	textureCoordinates->SetNumberOfComponents(2);
	textureCoordinates->SetName("TextureCoordinates");
	size_t nPoints = points->GetNumberOfPoints();
/*  // adaptation for netcdf-curtain for TERENO Demo
	double dist(0.0);
	for (size_t i = 0; i < nPoints; i++)
	{
		double coords[3];
		if ((i==0) || (i==173))
		{
			if (i==0) dist=0;
		}
		else
		{
			points->GetPoint(i-1, coords);
			GeoLib::Point* pnt = new GeoLib::Point(coords);
			points->GetPoint(i, coords);
			GeoLib::Point* pnt2 = new GeoLib::Point(coords);
			if (i<173)
				dist += sqrt(MathLib::sqrDist(pnt, pnt2));
			else
				dist -= sqrt(MathLib::sqrDist(pnt, pnt2));
		}
		points->GetPoint(i, coords);
		double x = MathLib::normalize(0, 8404, dist);
		double z = MathLib::normalize(-79.5, 1.5, coords[2]);
		float newcoords[2] = {x, z};
		textureCoordinates->InsertNextTuple(newcoords);
	}
*/
	for (size_t i = 0; i < nPoints; i++)
	{
		double coords[3];
		points->GetPoint(i, coords);
		float newcoords[2] = { MathLib::normalize(min.first, max.first, coords[0]),
		                       MathLib::normalize(min.second,max.second, coords[1])};
		textureCoordinates->InsertNextTuple(newcoords);
	}

	// put it all together
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	output->CopyStructure(input);
	output->GetPointData()->PassData(input->GetPointData());
	output->GetCellData()->PassData(input->GetCellData());
	output->GetPointData()->SetTCoords(textureCoordinates);

	return 1;
}

void VtkTextureOnSurfaceFilter::SetRaster(vtkImageAlgorithm* img,
                                          double x0, double y0,
                                          double scalingFactor)
{
	double range[2];
	img->Update();
	img->GetOutput()->GetPointData()->GetScalars()->GetRange(range);
	vtkSmartPointer<vtkImageShiftScale> scale = vtkSmartPointer<vtkImageShiftScale>::New();
	scale->SetInputConnection(img->GetOutputPort());
	scale->SetShift(-range[0]);
	scale->SetScale(255.0/(range[1]-range[0]));
	scale->SetOutputScalarTypeToUnsignedChar(); // Comment this out to get colored grayscale textures
	scale->Update();

	vtkTexture* texture = vtkTexture::New();
	texture->InterpolateOff();
	texture->RepeatOff();
	// texture->EdgeClampOn(); // does not work
	texture->SetInput(scale->GetOutput());
	this->SetTexture(texture);

	_origin = std::pair<float, float>(static_cast<float>(x0), static_cast<float>(y0));
	_scalingFactor = scalingFactor;
}

void VtkTextureOnSurfaceFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);
}
