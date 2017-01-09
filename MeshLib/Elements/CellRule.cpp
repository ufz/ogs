/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CellRule.h"

#include "MathLib/Vector3.h"
#include "MeshLib/Node.h"
#include "Element.h"
#include "FaceRule.h"

namespace MeshLib {

bool CellRule::testElementNodeOrder(const Element* e)
{
    const MathLib::Vector3 c (e->getCenterOfGravity());
    const unsigned nFaces (e->getNumberOfFaces());
    for (unsigned j=0; j<nFaces; ++j)
    {
        MeshLib::Element const*const face (e->getFace(j));
        // Node 1 is checked below because that way all nodes are used for the test
        // at some point, while for node 0 at least one node in every element
        // type would be used for checking twice and one wouldn't be checked at
        // all. (based on the definition of the _face_nodes variable)
        const MeshLib::Node x (*(face->getNode(1)));
        const MathLib::Vector3 cx (c, x);
        const double s = MathLib::scalarProduct(FaceRule::getSurfaceNormal(face), cx);
        delete face;
        if (s >= 0)
            return false;
    }
    return true;
}

} /* namespace */
