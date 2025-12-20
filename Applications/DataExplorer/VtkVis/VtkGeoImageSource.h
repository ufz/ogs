// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vtkSimpleImageToImageFilter.h>

#include <optional>

#include "GeoLib/Raster.h"
#include "VtkAlgorithmProperties.h"

namespace GeoLib
{
class Raster;
}

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

    // @brief Converts the source object into a Raster object
    static std::optional<GeoLib::Raster> convertToRaster(
        VtkGeoImageSource* const source);

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
