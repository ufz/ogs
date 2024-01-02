/**
 * \file
 * \author Lars Bilke
 * \date   2010-09-28
 * \brief  Implementation of the VtkGeoImageSource class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkGeoImageSource.h"

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

std::optional<GeoLib::Raster> VtkGeoImageSource::convertToRaster(
    VtkGeoImageSource* const source)
{
    int dims[3];
    source->GetOutput()->GetDimensions(dims);
    double origin[3];
    source->GetOutput()->GetOrigin(origin);
    double spacing[3];
    source->GetOutput()->GetSpacing(spacing);
    MathLib::Point3d const origin_pnt(
        std::array<double, 3>{{origin[0] - 0.5 * spacing[0],
                               origin[1] - 0.5 * spacing[1], origin[2]}});
    GeoLib::RasterHeader const header = {static_cast<std::size_t>(dims[0]),
                                         static_cast<std::size_t>(dims[1]),
                                         static_cast<std::size_t>(dims[2]),
                                         origin_pnt,
                                         spacing[0],
                                         -9999};

    vtkSmartPointer<vtkDataArray> const pixelData =
        vtkSmartPointer<vtkDataArray>(
            source->GetOutput()->GetPointData()->GetScalars());
    int const nTuple = pixelData->GetNumberOfComponents();
    if (nTuple < 1 || nTuple > 4)
    {
        ERR("VtkMeshConverter::convertImgToMesh(): Unsupported pixel "
            "composition!");
        return {};
    }

    std::vector<double> pix(header.n_cols * header.n_rows * header.n_depth, 0);
    for (std::size_t k = 0; k < header.n_depth; k++)
    {
        std::size_t const layer_idx = (k * header.n_rows * header.n_cols);
        for (std::size_t i = 0; i < header.n_rows; i++)
        {
            std::size_t const idx = i * header.n_cols + layer_idx;
            for (std::size_t j = 0; j < header.n_cols; j++)
            {
                double const* const colour = pixelData->GetTuple(idx + j);
                bool const visible = (nTuple == 2 || nTuple == 4)
                                         ? (colour[nTuple - 1] != 0)
                                         : true;
                if (!visible)
                {
                    pix[idx + j] = header.no_data;
                }
                else
                {
                    pix[idx + j] = (nTuple < 3) ? colour[0] :  // grey (+ alpha)
                                       (0.3 * colour[0] + 0.6 * colour[1] +
                                        0.1 * colour[2]);  // rgb(a)
                }
            }
        }
    }

    return std::make_optional<GeoLib::Raster>(header, pix.begin(), pix.end());
}
