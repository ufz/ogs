/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-21
 * \brief  Definition of the VtkCompositeColormapToImageFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"

/// @brief Applies a user adjustable color map to an image.
class VtkCompositeColormapToImageFilter : public VtkCompositeFilter
{
public:
    VtkCompositeColormapToImageFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeColormapToImageFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
};
