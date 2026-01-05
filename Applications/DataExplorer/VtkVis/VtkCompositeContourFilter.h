// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

/// @brief Visualisation of contour-lines/-planes within dense scalar fields.
/// In init() the threshold is first set to double min / max values. Set the
/// threshold later on via SetUserVectorProperty() to the actual data range.
class VtkCompositeContourFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeContourFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeContourFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
};
