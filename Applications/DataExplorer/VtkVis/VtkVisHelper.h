/**
 * \file
 * \author Lars Bilke
 * \date   2010-09-22
 * \brief  Definition of the VtkVisHelper class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKVISHELPER_H
#define VTKVISHELPER_H

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

#endif // VTKVISHELPER_H
