/**
 * \file
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Definition of the SizeDifferenceMetric class.
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
 * Calculates the quality of mesh elements based on the difference of element
 * size in comparison to the size of its neighbors.
 */
class SizeDifferenceMetric final : public ElementQualityMetric
{
public:
    explicit SizeDifferenceMetric(Mesh const& mesh);
    ~SizeDifferenceMetric() override = default;

    void calculateQuality() override;
};
}  // namespace MeshLib
