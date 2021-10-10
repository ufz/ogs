/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FaceRule.h"

#include "Element.h"
#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
bool FaceRule::testElementNodeOrder(Element const& e)
{
    return getSurfaceNormal(e)[2] < 0;
}

Eigen::Vector3d FaceRule::getFirstSurfaceVector(Element const& e)
{
    auto const a = Eigen::Map<Eigen::Vector3d const>(e.getNode(0)->getCoords());
    auto const b = Eigen::Map<Eigen::Vector3d const>(e.getNode(1)->getCoords());
    Eigen::Vector3d const v = a - b;
    return v;
}

Eigen::Vector3d FaceRule::getSecondSurfaceVector(Element const& e)
{
    auto const a = Eigen::Map<Eigen::Vector3d const>(e.getNode(1)->getCoords());
    auto const b = Eigen::Map<Eigen::Vector3d const>(e.getNode(2)->getCoords());
    Eigen::Vector3d const v = b - a;
    return v;
}

Eigen::Vector3d FaceRule::getSurfaceNormal(Element const& e)
{
    Eigen::Vector3d const u = getFirstSurfaceVector(e);
    Eigen::Vector3d const v = getSecondSurfaceVector(e);
    return u.cross(v);
}

}  // namespace MeshLib
