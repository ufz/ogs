// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief This filter colors the input by the points z-value.
class VtkCompositeColorByHeightFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeColorByHeightFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeColorByHeightFilter() override = default;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

protected:
};
