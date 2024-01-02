/**
 * \file
 * \author Karsten Rink
 * \date   2014-09-02
 * \brief  Definition of the RadiusEdgeRatioMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ElementQualityMetric.h"

namespace MeshToolsLib
{

/**
 * Calculates the quality of mesh elements based on the ratio between
 * radius of the smallest enclosing sphere and the shortest element edge
 */
struct RadiusEdgeRatioMetric final : public ElementQualityMetric
{
    using ElementQualityMetric::ElementQualityMetric;

    void calculateQuality() override;
};
}  // namespace MeshToolsLib
