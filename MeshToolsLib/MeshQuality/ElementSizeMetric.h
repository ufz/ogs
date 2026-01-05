// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ElementQualityMetric.h"

namespace MeshToolsLib
{

/**
 * Calculates the quality of mesh elements based on length/area/volume
 */
class ElementSizeMetric final : public ElementQualityMetric
{
public:
    using ElementQualityMetric::ElementQualityMetric;

    void calculateQuality() override;

private:
    std::size_t calc1dQuality();
    std::size_t calc2dOr3dQuality();
};
}  // namespace MeshToolsLib
