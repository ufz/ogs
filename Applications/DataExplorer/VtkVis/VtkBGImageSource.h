/**
 * \file
 * \author Karsten Rink
 * \date   2010-04-30
 * \brief  Definition of the VtkBGImageSource class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
