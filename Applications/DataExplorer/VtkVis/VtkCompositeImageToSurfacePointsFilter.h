// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

class VtkImageDataToSurfacePointsFilter;

class VtkCompositeImageToSurfacePointsFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeImageToSurfacePointsFilter(
        vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeImageToSurfacePointsFilter() override = default;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;
};
