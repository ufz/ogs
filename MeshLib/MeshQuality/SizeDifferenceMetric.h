/**
 * \file   SizeDifferenceMetric.h
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Definition of the SizeDifferenceMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
class SizeDifferenceMetric : public ElementQualityMetric
{
public:
    SizeDifferenceMetric(Mesh const& mesh);
    virtual ~SizeDifferenceMetric() = default;

    virtual void calculateQuality ();
};
}
