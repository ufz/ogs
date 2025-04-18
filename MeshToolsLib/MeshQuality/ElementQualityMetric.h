/**
 * \file
 * \author Thomas Fischer
 * \date   2010-12-08
 * \brief  Definition of the ElementQualityMetricBase class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Histogram.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace MeshToolsLib
{

/**
 * Base class for calculating the quality of mesh element based on a given
 * metric
 */
class ElementQualityMetric
{
public:
    explicit ElementQualityMetric(MeshLib::Mesh const& mesh);

    virtual ~ElementQualityMetric() = default;

    /// Calculates the quality metric for each element of the mesh
    virtual void calculateQuality() = 0;

    /// Returns the result vector
    std::vector<double> const& getElementQuality() const;

    /// Returns a histogram of the quality vector separated into the given
    /// number of bins. If no number of bins is specified, one will be
    /// calculated based on the Sturges criterium.
    virtual BaseLib::Histogram<double> getHistogram(
        std::size_t n_bins = 0) const;

protected:
    double _min = std::numeric_limits<double>::max();
    double _max = 0;
    MeshLib::Mesh const& _mesh;
    std::vector<double> _element_quality_metric;
};
}  // namespace MeshToolsLib
