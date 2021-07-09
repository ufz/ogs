/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Definition of the AngleSkewMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ElementQualityMetric.h"

namespace MeshLib
{

/**
 * Calculates the quality of mesh elements based on the EquiAngleSkew measure
 */
struct AngleSkewMetric final : public ElementQualityMetric
{
    using ElementQualityMetric::ElementQualityMetric;

    void calculateQuality() override;
};
}  // namespace MeshLib
