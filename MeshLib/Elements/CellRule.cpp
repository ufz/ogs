/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CellRule.h"

#include "Element.h"
#include "FaceRule.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
bool CellRule::testElementNodeOrder(Element const& e)
{
    Eigen::Vector3d const cc =
        Eigen::Map<Eigen::Vector3d const>(getCenterOfGravity(e).getCoords());
    const unsigned nFaces(e.getNumberOfFaces());
    for (unsigned j = 0; j < nFaces; ++j)
    {
        MeshLib::Element const* const face(e.getFace(j));
        // Node 1 is checked below because that way all nodes are used for the
        // test at some point, while for node 0 at least one node in every
        // element type would be used for checking twice and one wouldn't be
        // checked at all. (based on the definition of the _face_nodes variable)
        auto const x =
            Eigen::Map<Eigen::Vector3d const>(face->getNode(1)->getCoords());
        Eigen::Vector3d const cx = x - cc;
        const double s = FaceRule::getSurfaceNormal(*face).dot(cx);
        delete face;
        if (s >= 0)
        {
            return false;
        }
    }
    return true;
}

}  // namespace MeshLib
