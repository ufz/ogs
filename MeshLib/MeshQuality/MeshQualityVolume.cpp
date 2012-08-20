/*
 * MeshQualityVolume.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#include "MeshQualityVolume.h"

#include <iostream>

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
		MshElemType::type elem_type (elem->getType());
		if (elem_type == MshElemType::EDGE
		    || elem_type == MshElemType::TRIANGLE
		    || elem_type == MshElemType::QUAD)
		{
            _mesh_quality_measure[k] = -1.0;
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

	std::cout << "MeshQualityVolume::check() minimum: " << _min
	          << ", max_volume: " << _max << std::endl;
	if (error_count > 0)
		std::cout << "Warning: " << error_count << " elements with zero volume found." <<
		std::endl;
}

} // end namespace MeshLib
