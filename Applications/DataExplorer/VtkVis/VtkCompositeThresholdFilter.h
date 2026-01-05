// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

/// @brief Visualises only parts of meshes that are above/below/within given thresholds.
/// In init() the threshold is first set to double min / max values. Set the
/// threshold later on via SetUserVectorProperty() to the actual data range.
class VtkCompositeThresholdFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeThresholdFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeThresholdFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
};
