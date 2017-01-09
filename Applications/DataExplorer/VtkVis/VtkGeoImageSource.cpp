/**
 * \file
 * \author Lars Bilke
 * \date   2010-09-28
 * \brief  Implementation of the VtkGeoImageSource class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkGeoImageSource.h"

//#include "OGSRaster.h"
#include "VtkRaster.h"

#include <vtkFloatArray.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageShiftScale.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>

vtkStandardNewMacro(VtkGeoImageSource);

// This function is copied from the vtkSimpleImageFilterExample.cpp
//
// The switch statement in Execute will call this method with
// the appropriate input type (IT). Note that this example assumes
// that the output data type is the same as the input data type.
// This is not always the case.
template <class IT>
void vtkSimpleImageFilterExampleExecute(vtkImageData* input,
                                        vtkImageData* output,
                                        IT* inPtr, IT* outPtr)
{
    int dims[3];
    input->GetDimensions(dims);
    if (input->GetScalarType() != output->GetScalarType())
    {
        vtkGenericWarningMacro(<< "Execute: input ScalarType, " << input->GetScalarType()
                               << ", must match out ScalarType " << output->GetScalarType());
        return;
    }
    // HACK LB Multiply by number of scalar components due to RGBA values ?????
    int size = dims[0] * dims[1] * dims[2] * input->GetNumberOfScalarComponents();

    for(int i = 0; i < size; i++)
        outPtr[i] = inPtr[i];
}

VtkGeoImageSource::VtkGeoImageSource()
: _imageSource(NULL), _x0(0), _y0(0), _z0(0), _spacing(1)
{
}

VtkGeoImageSource::~VtkGeoImageSource()
{
    if(_imageSource) _imageSource->Delete();
}

void VtkGeoImageSource::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

bool VtkGeoImageSource::readImage(const QString &filename)
{
    vtkImageAlgorithm* img (VtkRaster::loadImage(filename.toStdString(), _x0, _y0, _spacing));
    if (img == nullptr)
        return false;
    this->setImage(img, filename, _x0, _y0, _spacing);
    return true;
}

void VtkGeoImageSource::setImage(vtkImageAlgorithm* image, const QString &name, double x0, double y0, double spacing)
{
    this->_imageSource = image;
    this->SetInputConnection(_imageSource->GetOutputPort());
    this->SetName(name);
    _x0 = x0; _y0 = y0; _z0 = -10; _spacing = spacing;

    this->GetOutput()->SetOrigin(_x0, _y0, _z0);
    this->GetOutput()->SetSpacing(_spacing, _spacing, _spacing);
}

vtkImageData* VtkGeoImageSource::getImageData()
{
    return this->_imageSource->GetImageDataInput(0);
}

void VtkGeoImageSource::getOrigin(double origin[3]) const
{
    origin[0] = this->_x0;
    origin[1] = this->_y0;
    origin[2] = this->_z0;
}

double VtkGeoImageSource::getSpacing() const
{
    return this->_spacing;
}

void VtkGeoImageSource::getRange(double range[2])
{
    this->_imageSource->Update();
    _imageSource->GetOutput()->GetPointData()->GetArray(0)->GetRange(range);
}

void VtkGeoImageSource::SimpleExecute(vtkImageData* input, vtkImageData* output)
{
    vtkDebugMacro(<< "Executing VtkGeoImageSource")
    void* inPtr = input->GetScalarPointer();
    void* outPtr = output->GetScalarPointer();
    switch(output->GetScalarType())
    {
        // This is simply a #define for a big case list.
        // It handles all data types that VTK supports.
        vtkTemplateMacro(vtkSimpleImageFilterExampleExecute(
                                 input, output, (VTK_TT*)(inPtr), (VTK_TT*)(outPtr)));
    default:
        vtkGenericWarningMacro("Execute: Unknown input ScalarType");
        return;
    }
}

void VtkGeoImageSource::SetUserProperty( QString name, QVariant value )
{
    Q_UNUSED(name);
    Q_UNUSED(value);
}
