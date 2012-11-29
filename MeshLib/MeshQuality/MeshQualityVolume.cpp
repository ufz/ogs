/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MeshQualityVolume.cpp
 *
 *  Created on 2011-03-03 by Thomas Fischer
 */

#include "MeshQualityVolume.h"

#include <iostream>
#include <limits>

namespace MeshLib
{

MeshQualityVolume::MeshQualityVolume(Mesh const* const mesh) :
	MeshQualityChecker(mesh)
{ }

void MeshQualityVolume::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::Element*>& elements(_mesh->getElements());

	size_t error_count(0);
	size_t nElements (_mesh->getNElements());

	for (size_t k(0); k < nElements; k++)
	{
		const Element* elem (elements[k]);
		MshElemType::type elem_type (elem->getGeomType());
		if (elem_type == MshElemType::EDGE
		    || elem_type == MshElemType::TRIANGLE
		    || elem_type == MshElemType::QUAD)
		{
            _mesh_quality_measure[k] = 0.0;
            continue;
        }

        double volume (elem->getContent());
        if (volume > _max)
            _max = volume;
        if (volume < sqrt(fabs(std::numeric_limits<double>::min()))) {
			errorMsg(elem, k);
			error_count++;
		} else if (volume < _min)
            _min = volume;
        _mesh_quality_measure[k] = volume;
	}

	INFO ("MeshQualityVolume::check() minimum: %f, max_volume: %f", _min, _max);
	if (error_count > 0)
		WARN ("Warning: %d elements with zero volume found.", error_count);
}

} // end namespace MeshLib
