/**
 * \file VtkBGImageSource.cpp
 * 30/04/2010 KR Initial implementation
 *
 */

// ** INCLUDES **
#include "VtkBGImageSource.h"
#include "VtkVisHelper.h"

#include <vtkImageAlgorithm.h>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTexture.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageShiftScale.h>

vtkStandardNewMacro(VtkBGImageSource);
vtkCxxRevisionMacro(VtkBGImageSource, "$Revision$");

VtkBGImageSource::VtkBGImageSource() : _origin(0,0), _cellsize(1)
{
}

VtkBGImageSource::~VtkBGImageSource()
{
}

void VtkBGImageSource::SetRaster(vtkImageAlgorithm *img, double x0, double y0, double scalingFactor)
{
	double range[2];
	img->Update();
	img->GetOutput()->GetPointData()->GetScalars()->GetRange(range);
	vtkSmartPointer<vtkImageShiftScale> scale = vtkSmartPointer<vtkImageShiftScale>::New();
	scale->SetInputConnection(img->GetOutputPort());
	scale->SetShift(-range[0]);
	scale->SetScale(255.0/(range[1]-range[0]));
	scale->SetOutputScalarTypeToUnsignedChar();
	scale->Update();

	vtkImageData* imageData = scale->GetOutput();
	int dims[3];
	imageData->GetDimensions(dims);
	vtkTexture* texture = vtkTexture::New();
	texture->InterpolateOff();
	texture->RepeatOff();
	// texture->EdgeClampOn(); // does not work
	texture->SetInput(imageData);
	this->SetTexture(texture);

	_origin = std::pair<float, float>(static_cast<float>(x0), static_cast<float>(y0));
	_cellsize = scalingFactor;


	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
	plane->SetOrigin( _origin.first, _origin.second, -1 );
	plane->SetPoint1( _origin.first + dims[0] * _cellsize, _origin.second, -1 );
	plane->SetPoint2( _origin.first,_origin.second + dims[1] * _cellsize, -1 );

	this->SetInputConnection(0, plane->GetOutputPort(0));
	this->SetTexture(texture);
}

void VtkBGImageSource::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);
}
