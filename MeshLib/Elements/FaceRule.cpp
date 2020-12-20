/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

Eigen::Vector3d FaceRule::getFirstSurfaceVector(Element const* const e)
{
    auto const a =
        Eigen::Map<Eigen::Vector3d const>(e->getNode(0)->getCoords());
    auto const b =
        Eigen::Map<Eigen::Vector3d const>(e->getNode(1)->getCoords());
    Eigen::Vector3d const v = a - b;
    return v;
}

Eigen::Vector3d FaceRule::getSecondSurfaceVector(Element const* const e)
{
    auto const a =
        Eigen::Map<Eigen::Vector3d const>(e->getNode(1)->getCoords());
    auto const b =
        Eigen::Map<Eigen::Vector3d const>(e->getNode(2)->getCoords());
    Eigen::Vector3d const v = b - a;
    return v;
}

MathLib::Vector3 FaceRule::getSurfaceNormal(const Element* e)
{
    Eigen::Vector3d const u = getFirstSurfaceVector(e);
    Eigen::Vector3d const v = getSecondSurfaceVector(e);
    Eigen::Vector3d const normal = u.cross(v);
    return MathLib::Vector3{normal[0], normal[1], normal[2]};
}

}  // namespace MeshLib
