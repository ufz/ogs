// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief Puts a texture on an object (and converts it into a vtkPolyData if necessary).
class VtkCompositeTextureOnSurfaceFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeTextureOnSurfaceFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeTextureOnSurfaceFilter() override;

    void init() override;

    void SetUserProperty(QString name, QVariant value) override;

private:
};
