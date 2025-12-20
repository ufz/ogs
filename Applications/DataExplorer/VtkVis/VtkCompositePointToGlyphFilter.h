// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief Converts point data to scalar-scaled spheres.
class VtkCompositePointToGlyphFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositePointToGlyphFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositePointToGlyphFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

private:
    vtkSphereSource* _glyphSource;
};
