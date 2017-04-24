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
    return getSurfaceNormal(e)[2] < 0;
}

MathLib::Vector3 FaceRule::getFirstSurfaceVector(Element const* const e)
{
    Node* const* const _nodes = e->getNodes();
    return {*_nodes[1], *_nodes[0]};
}

MathLib::Vector3 FaceRule::getSecondSurfaceVector(Element const* const e)
{
    Node* const* const _nodes = e->getNodes();
    return {*_nodes[1], *_nodes[2]};
}

MathLib::Vector3 FaceRule::getSurfaceNormal(const Element* e)
{
    const MathLib::Vector3 u = getFirstSurfaceVector(e);
    const MathLib::Vector3 v = getSecondSurfaceVector(e);
    return MathLib::crossProduct(u, v);
}

} /* namespace */

