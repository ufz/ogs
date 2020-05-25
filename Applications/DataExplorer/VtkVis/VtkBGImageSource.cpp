/**
 * \file
 * \author Karsten Rink
 * \date   2010-04-30
 * \brief  Implementation of the VtkBGImageSource class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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

VtkBGImageSource::VtkBGImageSource() = default;
VtkBGImageSource::~VtkBGImageSource() = default;

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
    texture->SetInputData(imageData);
    this->SetTexture(texture);

    origin_ = std::pair<float, float>(static_cast<float>(x0), static_cast<float>(y0));
    cellsize_ = scalingFactor;


    vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
    plane->SetOrigin( origin_.first, origin_.second, -1 );
    plane->SetPoint1( origin_.first + dims[0] * cellsize_, origin_.second, -1 );
    plane->SetPoint2( origin_.first,origin_.second + dims[1] * cellsize_, -1 );

    this->SetInputConnection(0, plane->GetOutputPort(0));
    this->SetTexture(texture);
}

void VtkBGImageSource::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);
}
