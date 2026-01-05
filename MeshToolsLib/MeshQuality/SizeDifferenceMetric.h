// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ElementQualityMetric.h"

namespace MeshToolsLib
{

/**
 * Calculates the quality of mesh elements based on the difference of element
 * size in comparison to the size of its neighbors.
 */
struct SizeDifferenceMetric final : public ElementQualityMetric
{
    using ElementQualityMetric::ElementQualityMetric;

    void calculateQuality() override;
};
}  // namespace MeshToolsLib
