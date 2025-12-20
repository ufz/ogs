// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"
#include "GeoLib/GeoType.h"

class vtkThreshold;

/// @brief Highlights a single GeoObject
class VtkCompositeGeoObjectFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeGeoObjectFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeGeoObjectFilter() override;

    void init() override;

    /// @brief Sets user properties.
    void SetUserProperty(QString name, QVariant value) override
    {
        Q_UNUSED(name);
        Q_UNUSED(value);
    }

    void SetIndex(std::size_t idx);

private:
    GeoLib::GEOTYPE _type;
    vtkThreshold* _threshold;
};
