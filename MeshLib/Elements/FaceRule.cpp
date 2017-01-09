/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FaceRule.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"
#include "Element.h"

namespace MeshLib
{

bool FaceRule::testElementNodeOrder(const Element* e)
{
    MathLib::Vector3 up_vec (0,0,1);
    return (MathLib::scalarProduct(getSurfaceNormal(e), up_vec) < 0) ? true : false;
}

MathLib::Vector3 FaceRule::getSurfaceNormal(const Element* e)
{
    Node * const * _nodes = e->getNodes();
    const MathLib::Vector3 u (*_nodes[1], *_nodes[0]);
    const MathLib::Vector3 v (*_nodes[1], *_nodes[2]);
    return MathLib::crossProduct(u,v);
}

} /* namespace */

