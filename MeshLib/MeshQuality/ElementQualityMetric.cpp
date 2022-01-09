/**
 * \file
 * \author Thomas Fischer
 * \date   2010-12-08
 * \brief  Implementation of the ElementQualityMetricBase class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementQualityMetric.h"

#include <cmath>

#include "MeshLib/Node.h"

namespace MeshLib
{
ElementQualityMetric::ElementQualityMetric(Mesh const& mesh) : _mesh(mesh)
{
    _element_quality_metric.resize(_mesh.getNumberOfElements(), -1.0);
}

BaseLib::Histogram<double> ElementQualityMetric::getHistogram(
    std::size_t n_bins) const
{
    if (n_bins == 0)
    {
        n_bins = static_cast<std::size_t>(
            1 +
            3.3 * std::log(static_cast<float>((_mesh.getNumberOfElements()))));
    }

    return BaseLib::Histogram<double>(getElementQuality(), n_bins, true);
}

std::vector<double> const& ElementQualityMetric::getElementQuality() const
{
    return _element_quality_metric;
}
}  // namespace MeshLib
