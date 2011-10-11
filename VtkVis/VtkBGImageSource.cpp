/**
 * \file VtkBGImageSource.cpp
 * 30/04/2010 KR Initial implementation
 *
 */

// ** INCLUDES **
#include "VtkBGImageSource.h"
#include "VtkVisHelper.h"

#include <QImage>

#include "vtkObjectFactory.h"
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTexture.h>

vtkStandardNewMacro(VtkBGImageSource);
vtkCxxRevisionMacro(VtkBGImageSource, "$Revision$");

VtkBGImageSource::VtkBGImageSource() : _origin(0,0), _cellsize(1)
{
}

VtkBGImageSource::~VtkBGImageSource()
{
}

void VtkBGImageSource::SetImage(vtkTexture* texture)
{
	this->SetTexture(texture);
}

void VtkBGImageSource::SetRaster(QImage &img)
{
	vtkTexture* texture = VtkVisHelper::QImageToVtkTexture(img);

	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
	plane->SetOrigin( _origin.first, _origin.second, -1 );
	plane->SetPoint1( _origin.first + img.width() * _cellsize, _origin.second, -1 );
	plane->SetPoint2( _origin.first,_origin.second + img.height() * _cellsize, -1 );

	this->SetInputConnection(0, plane->GetOutputPort(0));
	this->SetTexture(texture);
}

void VtkBGImageSource::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);
}
