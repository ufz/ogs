/**
 * \file
 * \author Lars Bilke
 * \date   2010-09-22
 * \brief  Implementation of the VtkVisHelper class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkVisHelper.h"

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTexture.h>
#include <vtkUnsignedCharArray.h>

#include <QImage>

vtkImageData* VtkVisHelper::QImageToVtkImageData(QImage &img)
{
    std::size_t imgWidth = img.width(), imgHeight = img.height();
    vtkSmartPointer<vtkUnsignedCharArray> data = vtkSmartPointer<vtkUnsignedCharArray>::New();
    data->SetNumberOfComponents(3);
    data->SetNumberOfTuples( imgWidth * imgHeight );

    for (std::size_t j = 0; j < imgHeight; j++)
        for (std::size_t i = 0; i < imgWidth; i++)
        {
            QRgb pix = img.pixel(i,j);
            const float color[3] = { static_cast<float>(qRed(pix)),
                                     static_cast<float>(qGreen(pix)),
                                     static_cast<float>(qBlue(pix))
                                   };
            data->SetTuple(j * imgWidth + i, color);
        }

    vtkImageData* imgData = vtkImageData::New();
    imgData->SetExtent(0, imgWidth - 1, 0, imgHeight - 1, 0, 0);
    imgData->SetOrigin(0, 0, 0);
    imgData->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    imgData->GetPointData()->SetScalars(data);

    return imgData;
}

vtkTexture* VtkVisHelper::QImageToVtkTexture(QImage &img)
{
    vtkSmartPointer<vtkImageData> imgData = QImageToVtkImageData(img);

    vtkTexture* texture = vtkTexture::New();
    texture->InterpolateOff();
    texture->RepeatOff();
    //texture->EdgeClampOff();
    //texture->SetBlendingMode(0);
    texture->SetInputData(imgData);

    return texture;
}
