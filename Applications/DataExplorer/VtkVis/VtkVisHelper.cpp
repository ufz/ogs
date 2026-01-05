// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "VtkVisHelper.h"

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTexture.h>
#include <vtkUnsignedCharArray.h>

#include <QImage>

vtkImageData* VtkVisHelper::QImageToVtkImageData(QImage& img)
{
    std::size_t imgWidth = img.width();
    std::size_t imgHeight = img.height();
    vtkSmartPointer<vtkUnsignedCharArray> data =
        vtkSmartPointer<vtkUnsignedCharArray>::New();
    data->SetNumberOfComponents(3);
    data->SetNumberOfTuples(imgWidth * imgHeight);

    for (std::size_t j = 0; j < imgHeight; j++)
    {
        for (std::size_t i = 0; i < imgWidth; i++)
        {
            QRgb pix = img.pixel(i, j);
            const float color[3] = {static_cast<float>(qRed(pix)),
                                    static_cast<float>(qGreen(pix)),
                                    static_cast<float>(qBlue(pix))};
            data->SetTuple(j * imgWidth + i, color);
        }
    }

    vtkImageData* imgData = vtkImageData::New();
    imgData->SetExtent(0, imgWidth - 1, 0, imgHeight - 1, 0, 0);
    imgData->SetOrigin(0, 0, 0);
    imgData->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    imgData->GetPointData()->SetScalars(data);

    return imgData;
}

vtkTexture* VtkVisHelper::QImageToVtkTexture(QImage& img)
{
    vtkSmartPointer<vtkImageData> imgData = QImageToVtkImageData(img);

    vtkTexture* texture = vtkTexture::New();
    texture->InterpolateOff();
    texture->RepeatOff();
    // texture->EdgeClampOff();
    // texture->SetBlendingMode(0);
    texture->SetInputData(imgData);

    return texture;
}
