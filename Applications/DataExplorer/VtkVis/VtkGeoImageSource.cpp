/**
 * \file
 * \author Lars Bilke
 * \date   2010-09-28
 * \brief  Implementation of the VtkGeoImageSource class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkGeoImageSource.h"

//#include "OGSRaster.h"
#include <vtkFloatArray.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageShiftScale.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

#include "VtkRaster.h"

vtkStandardNewMacro(VtkGeoImageSource);

// This function is copied from the vtkSimpleImageFilterExample.cpp
//
// The switch statement in Execute will call this method with
// the appropriate input type (IT). Note that this example assumes
// that the output data type is the same as the input data type.
// This is not always the case.
template <class IT>
void vtkSimpleImageFilterExampleExecute(vtkImageData* input,
                                        vtkImageData* output, IT* inPtr,
                                        IT* outPtr)
{
    int dims[3];
    input->GetDimensions(dims);
    if (input->GetScalarType() != output->GetScalarType())
    {
        vtkGenericWarningMacro(
            << "Execute: input ScalarType, " << input->GetScalarType()
            << ", must match out ScalarType " << output->GetScalarType());
        return;
    }
    // HACK LB Multiply by number of scalar components due to RGBA values ?????
    int size =
        dims[0] * dims[1] * dims[2] * input->GetNumberOfScalarComponents();

    for (int i = 0; i < size; i++)
    {
        outPtr[i] = inPtr[i];
    }
}

VtkGeoImageSource::VtkGeoImageSource() = default;

VtkGeoImageSource::~VtkGeoImageSource()
{
    if (_imageSource)
    {
        _imageSource->Delete();
    }
}

void VtkGeoImageSource::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

bool VtkGeoImageSource::readImage(const QString& filename)
{
    vtkImageAlgorithm* img(VtkRaster::loadImage(filename.toStdString()));
    if (img == nullptr)
    {
        return false;
    }
    this->setImage(img, filename);
    return true;
}

void VtkGeoImageSource::setImage(vtkImageAlgorithm* image, const QString& name)
{
    this->_imageSource = image;
    this->SetInputConnection(_imageSource->GetOutputPort());
    this->SetName(name);
}

vtkImageData* VtkGeoImageSource::getImageData()
{
    return this->_imageSource->GetImageDataInput(0);
}

void VtkGeoImageSource::SimpleExecute(vtkImageData* input, vtkImageData* output)
{
    vtkDebugMacro(<< "Executing VtkGeoImageSource");
    void* inPtr = input->GetScalarPointer();
    void* outPtr = output->GetScalarPointer();
    switch (output->GetScalarType())
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

void VtkGeoImageSource::SetUserProperty(QString name, QVariant value)
{
    Q_UNUSED(name);
    Q_UNUSED(value);
}
