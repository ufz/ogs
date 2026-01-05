// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

/// @brief Converts lines to tube-objects.
class VtkCompositeLineToTubeFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeLineToTubeFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeLineToTubeFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

private:

};
