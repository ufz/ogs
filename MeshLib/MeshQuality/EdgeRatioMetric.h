/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Definition of the EdgeRatioMetric class.
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
#include "MathLib/Point3d.h"

namespace MeshLib
{

/**
 * Calculates the quality of mesh elements based on the ratio between shortest and longest edge of an element
 */
class EdgeRatioMetric final : public ElementQualityMetric
{
public:
    explicit EdgeRatioMetric(Mesh const& mesh);
    ~EdgeRatioMetric() override = default;

    void calculateQuality() override;
};
}  // namespace MeshLib
