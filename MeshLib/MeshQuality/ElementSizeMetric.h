/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-17
 * \brief  Implementation of the AreaMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
 * Calculates the quality of mesh elements based on length/area/volume
 */
class ElementSizeMetric final : public ElementQualityMetric
{
public:
    using ElementQualityMetric::ElementQualityMetric;

    void calculateQuality() override;

private:
    std::size_t calc1dQuality();
    std::size_t calc2dQuality();
    std::size_t calc3dQuality();
};
}  // namespace MeshLib
