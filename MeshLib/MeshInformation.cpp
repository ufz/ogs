/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Implementation of the MeshInformation class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshInformation.h"
#include "Mesh.h"
#include "Elements/Element.h"

namespace MeshLib
{

const std::pair<int, int> MeshInformation::getValueBounds(const MeshLib::Mesh &mesh)
{
	boost::optional<MeshLib::PropertyVector<int> const&> materialIds
		= mesh.getProperties().getPropertyVector<int>("MaterialIDs");
	if (!materialIds) {
		INFO("Mesh does not contain a property \"MaterialIDs\".");
		return {std::numeric_limits<int>::max(),
		        std::numeric_limits<int>::max()};
	}
	if (materialIds->empty()) {
		INFO("Mesh does not contain values for the property \"MaterialIDs\".");
		return {std::numeric_limits<int>::max(),
		        std::numeric_limits<int>::max()};
	}
	auto mat_bounds = std::minmax_element(materialIds->cbegin(), materialIds->cend());
	return {*(mat_bounds.first), *(mat_bounds.second)};
}

const GeoLib::AABB<MeshLib::Node> MeshInformation::getBoundingBox(const MeshLib::Mesh &mesh)
{
	const std::vector<MeshLib::Node*> &nodes (mesh.getNodes());
	return GeoLib::AABB<MeshLib::Node>(nodes.begin(), nodes.end());
}

const std::array<unsigned, 7> MeshInformation::getNumberOfElementTypes(const MeshLib::Mesh &mesh)
{
	std::array<unsigned, 7> n_element_types = {{0, 0, 0, 0, 0, 0, 0}};
	const std::vector<MeshLib::Element*> &elements (mesh.getElements());
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

} //end MeshLib
