// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// ** INCLUDES **
#include <vtkTextureMapToPlane.h>

#include "VtkAlgorithmProperties.h"

class vtkImageAlgorithm;

/**
 * \brief Uses an image source to create a plane in the 3D with the given
 * image texture mapped on it.
 */
class VtkBGImageSource : public vtkTextureMapToPlane, public VtkAlgorithmProperties
{
public:
    /// Create new objects with New() because of VTKs object reference counting.
    static VtkBGImageSource* New();

    vtkTypeMacro(VtkBGImageSource, vtkTextureMapToPlane);

    /// Sets the raster/image to be used as a texture map
    void SetRaster(vtkImageAlgorithm *img, double x0, double y0, double scalingFactor);

    void SetUserProperty(QString name, QVariant value) override;

protected:
    VtkBGImageSource();
    ~VtkBGImageSource() override;

private:
    std::pair<double, double> _origin{0, 0};
    double _cellsize{1};
};
