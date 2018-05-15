/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"

class VtkImageDataToSurfacePointsFilter;

class VtkCompositeImageToSurfacePointsFilter : public VtkCompositeFilter
{
public:
    VtkCompositeImageToSurfacePointsFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeImageToSurfacePointsFilter() = default;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;
};
