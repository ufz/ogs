/**
 * \file   ElementQualityMetric.cpp
 * \author Thomas Fischer
 * \date   2010-12-08
 * \brief  Implementation of the ElementQualityMetricBase class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementQualityMetric.h"
#include "Node.h"
#include "Point.h"
#include <cmath>
#include <iostream>

namespace MeshLib
{
ElementQualityMetric::ElementQualityMetric(Mesh const* const mesh) :
	_min (std::numeric_limits<double>::max()), _max (std::numeric_limits<double>::min()), _mesh (mesh)
{
	if (_mesh)
		_element_quality_metric.resize (_mesh->getNElements(), -1.0);
}

BaseLib::Histogram<double> ElementQualityMetric::getHistogram (size_t nclasses) const
{
	if (nclasses == 0) {
		// simple suggestion: number of classes with Sturges criterion
		nclasses = static_cast<size_t>(1 + 3.3 * log (static_cast<float>((_mesh->getNElements()))));
	}

	return BaseLib::Histogram<double>(getElementQuality(), nclasses, true);
}

void ElementQualityMetric::errorMsg (const Element* elem, size_t idx) const
{
	ERR ("Error in MeshQualityChecker::check() - Calculated value of element is below double precision minimum.");
	ERR ("Points of %s-Element %d: ", MeshElemType2String(elem->getGeomType()).c_str(), idx);
	for (size_t i(0); i < elem->getNNodes(); i++)
	{
		const double* coords = elem->getNode(i)->getCoords();
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
