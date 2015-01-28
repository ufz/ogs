/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CellRule.h"

#include "logog/include/logog.hpp"

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Node.h"
#include "Face.h"

namespace MeshLib
{

const unsigned CellRule::dimension;

bool CellRule::testElementNodeOrder(const Element* e)
{
	const MathLib::Vector3 c (e->getCenterOfGravity());
	const unsigned nFaces (e->getNFaces());
	for (unsigned j=0; j<nFaces; ++j)
	{
		MeshLib::Face const*const face (dynamic_cast<const MeshLib::Face*>(e->getFace(j)));
		// Node 1 is checked below because that way all nodes are used for the test
		// at some point, while for node 0 at least one node in every element
		// type would be used for checking twice and one wouldn't be checked at
		// all. (based on the definition of the _face_nodes variable)
		const MeshLib::Node x (*(face->getNode(1)));
		const MathLib::Vector3 cx (c, x);
		const double s = MathLib::scalarProduct(face->getSurfaceNormal(), cx);
		delete face;
		if (s >= 0)
			return false;
	}
	return true;
}

} // end namespace MeshLib
