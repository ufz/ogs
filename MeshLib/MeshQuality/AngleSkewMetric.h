/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Definition of the AngleSkewMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
class AngleSkewMetric final : public ElementQualityMetric
{
public:
    explicit AngleSkewMetric(Mesh const& mesh);

    void calculateQuality() override;
};
}  // namespace MeshLib
