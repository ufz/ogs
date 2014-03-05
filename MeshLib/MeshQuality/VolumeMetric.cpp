/**
 * \file   VolumeMetric.cpp
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Implementation of the VolumeMetric class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VolumeMetric.h"

#include <iostream>
#include <limits>

namespace MeshLib
{

VolumeMetric::VolumeMetric(Mesh const* const mesh) :
	ElementQualityMetric(mesh)
{ }

void VolumeMetric::calculateQuality()
{
	// get all elements of mesh
	const std::vector<MeshLib::Element*>& elements(_mesh->getElements());

	size_t error_count(0);
	size_t nElements (_mesh->getNElements());

	for (size_t k(0); k < nElements; k++)
	{
		const Element* elem (elements[k]);
		if (elem->getDimension()<3)
		{
            _element_quality_metric[k] = 0.0;
            continue;
        }

        double volume (elem->getContent());
        if (volume > _max)
            _max = volume;
        if (volume < sqrt(fabs(std::numeric_limits<double>::epsilon()))) {
			errorMsg(elem, k);
			error_count++;
		} else if (volume < _min)
            _min = volume;
        _element_quality_metric[k] = volume;
	}

	INFO ("MeshQualityVolume::check() minimum: %f, max_volume: %f", _min, _max);
	if (error_count > 0)
		WARN ("Warning: %d elements with zero volume found.", error_count);
}

} // end namespace MeshLib
