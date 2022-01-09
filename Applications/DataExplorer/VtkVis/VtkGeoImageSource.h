/**
 * \file
 * \author Lars Bilke
 * \date   2010-09-28
 * \brief  Definition of the VtkGeoImageSource class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkAlgorithmProperties.h"
#include <vtkSimpleImageToImageFilter.h>

class QString;
class QPointF;
class QImage;
class vtkQImageToImageSource;
class vtkImageShiftScale;
class vtkImageData;

/**
 * \brief The VtkVisPipeline source object of a geo-referenced image (file).
 */
class VtkGeoImageSource : public vtkSimpleImageToImageFilter, public VtkAlgorithmProperties
{
public:
    /// @brief Create new objects with New() because of VTKs reference counting.
    static VtkGeoImageSource* New();

    vtkTypeMacro(VtkGeoImageSource, vtkSimpleImageToImageFilter);

    /// @brief Prints information about itself.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    VtkGeoImageSource(const VtkGeoImageSource&) = delete;
    void operator=(const VtkGeoImageSource&) = delete;

    /// @brief Returns the ImageData object.
    vtkImageData* getImageData();

    /// @brief Reads an image from file.
    bool readImage(const QString &filename);

    /// @brief Imports an existing image object.
    void setImage(vtkImageAlgorithm* image, const QString& name);

    void SetUserProperty(QString name, QVariant value) override;

protected:
    /// @brief Constructor.
    VtkGeoImageSource();

    /// @brief Destructor.
    ~VtkGeoImageSource() override;

    /// @brief Filter execution.
    void SimpleExecute(vtkImageData* input, vtkImageData* output) override;

private:
    vtkImageAlgorithm* _imageSource{nullptr};
};
