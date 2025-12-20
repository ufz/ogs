// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FaceRule.h"

#include <Eigen/Geometry>

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
    Eigen::Vector3d const v =
        e.getNode(0)->asEigenVector3d() - e.getNode(1)->asEigenVector3d();
    return v;
}

Eigen::Vector3d FaceRule::getSecondSurfaceVector(Element const& e)
{
    Eigen::Vector3d const v =
        e.getNode(2)->asEigenVector3d() - e.getNode(1)->asEigenVector3d();
    return v;
}

Eigen::Vector3d FaceRule::getSurfaceNormal(Element const& e)
{
    Eigen::Vector3d const u = getFirstSurfaceVector(e);
    Eigen::Vector3d const v = getSecondSurfaceVector(e);
    return u.cross(v);
}

}  // namespace MeshLib
