// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ElementQualityMetric.h"
#include "MathLib/Point3d.h"

namespace MeshToolsLib
{

/**
 * Calculates the quality of mesh elements based on the ratio between shortest
 * and longest edge of an element
 */
struct EdgeRatioMetric final : public ElementQualityMetric
{
    using ElementQualityMetric::ElementQualityMetric;

    void calculateQuality() override;
};
}  // namespace MeshToolsLib
