/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Implementation of the MeshInformation class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshInformation.h"
#include "Mesh.h"
#include "Elements/Element.h"


const std::pair<unsigned, unsigned> MeshInformation::getValueBounds(MeshLib::Mesh const*const mesh)
{
	const std::vector<MeshLib::Element*> elements (mesh->getElements());
	const auto minmax = std::minmax_element(elements.cbegin(), elements.cend(),
        [](MeshLib::Element const*const a, MeshLib::Element const*const b)
            {
                return a->getValue() < b->getValue();
        });
	return std::make_pair<unsigned, unsigned>((*minmax.first)->getValue(), (*minmax.second)->getValue());
}

const GeoLib::AABB<MeshLib::Node> MeshInformation::getBoundingBox(MeshLib::Mesh const*const mesh)
{
	const std::vector<MeshLib::Node*> nodes (mesh->getNodes());
	return GeoLib::AABB<MeshLib::Node>(nodes.begin(), nodes.end());
}

const std::array<unsigned, 7> MeshInformation::getNumberOfElementTypes(MeshLib::Mesh const*const mesh)
{
	std::array<unsigned, 7> n_element_types = { 0, 0, 0, 0, 0, 0, 0};
	const std::vector<MeshLib::Element*> elements (mesh->getElements());
	for (auto it = elements.begin(); it != elements.end(); ++it)
	{
		MeshElemType t = (*it)->getGeomType();
		if (t == MeshElemType::LINE) n_element_types[0]++;
		if (t == MeshElemType::TRIANGLE) n_element_types[1]++;
		if (t == MeshElemType::QUAD) n_element_types[2]++;
		if (t == MeshElemType::TETRAHEDRON) n_element_types[3]++;
		if (t == MeshElemType::HEXAHEDRON) n_element_types[4]++;
		if (t == MeshElemType::PYRAMID) n_element_types[5]++;
		if (t == MeshElemType::PRISM) n_element_types[6]++;
	}
	return n_element_types;
}