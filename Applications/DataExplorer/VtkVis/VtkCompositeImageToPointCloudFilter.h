// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

class VtkImageDataToPointCloudFilter;

class VtkCompositeImageToPointCloudFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeImageToPointCloudFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeImageToPointCloudFilter() override = default;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;
};
