/**
 * \file
 * \author Thomas Fischer
 * \date   2010-12-08
 * \brief  Implementation of the ElementQualityMetricBase class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
ElementQualityMetric::ElementQualityMetric(Mesh const& mesh) :
    min_ (std::numeric_limits<double>::max()), max_ (0), mesh_ (mesh)
{
    element_quality_metric_.resize (mesh_.getNumberOfElements(), -1.0);
}

BaseLib::Histogram<double> ElementQualityMetric::getHistogram (std::size_t n_bins) const
{
    if (n_bins == 0)
    {
        n_bins = static_cast<std::size_t>(
            1 + 3.3 * log(static_cast<float>((mesh_.getNumberOfElements()))));
    }

    return BaseLib::Histogram<double>(getElementQuality(), n_bins, true);
}

std::vector<double> const& ElementQualityMetric::getElementQuality () const
{
    return element_quality_metric_;
}
}  // namespace MeshLib
