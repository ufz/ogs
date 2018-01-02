/**
 * \file   ElementQualityMetric.cpp
 * \author Thomas Fischer
 * \date   2010-12-08
 * \brief  Implementation of the ElementQualityMetricBase class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    _min (std::numeric_limits<double>::max()), _max (0), _mesh (mesh)
{
    _element_quality_metric.resize (_mesh.getNumberOfElements(), -1.0);
}

BaseLib::Histogram<double> ElementQualityMetric::getHistogram (std::size_t n_bins) const
{
    if (n_bins == 0)
        n_bins = static_cast<std::size_t>(1 + 3.3 * log (static_cast<float>((_mesh.getNumberOfElements()))));

    return BaseLib::Histogram<double>(getElementQuality(), n_bins, true);
}

void ElementQualityMetric::errorMsg (Element const& elem, std::size_t idx) const
{
    ERR ("Error in MeshQualityChecker::check() - Calculated value of element is below double precision minimum.");
    ERR ("Points of %s-Element %d: ", MeshElemType2String(elem.getGeomType()).c_str(), idx);
    for (std::size_t i(0); i < elem.getNumberOfBaseNodes(); i++)
    {
        const double* coords = elem.getNode(i)->getCoords();
        ERR ("\t Node %d: (%f, %f, %f)", i, coords[0], coords[1], coords[2]);
    }
}

std::vector<double> const& ElementQualityMetric::getElementQuality () const
{
    return _element_quality_metric;
}

double ElementQualityMetric::getMinValue() const
{
    return _min;
}

double ElementQualityMetric::getMaxValue() const
{
    return _max;
}

} // end namespace MeshLib
