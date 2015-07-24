/**
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "MeshNodesOnPoint.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{

MeshNodesOnPoint::MeshNodesOnPoint(MeshLib::Mesh const& mesh,
        GeoLib::Point const& pnt, double epsilon_radius, bool search_all_nodes)
: _mesh(mesh), _pnt(pnt)
{
	for (auto node : mesh.getNodes())
	{
		double len(sqrt(MathLib::sqrDist(pnt, *node)));
		if (len < epsilon_radius)
			_msh_node_ids.push_back(node->getID());
	}

}

} // end namespace MeshGeoToolsLib

