// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

class QImage;
class vtkTexture;
class vtkImageData;

/**
 * \brief Some data conversion helper functions.
 */
class VtkVisHelper
{
public:
    /// @brief Converts a QImage to vtkImageData.
    static vtkImageData* QImageToVtkImageData(QImage &img);

    /// @brief Converts a QImage-object into a vtkTexture-object
    static vtkTexture* QImageToVtkTexture(QImage &img);
};
