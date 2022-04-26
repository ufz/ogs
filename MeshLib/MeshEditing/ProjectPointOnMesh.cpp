/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoLib/AnalyticalGeometry.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"

namespace MeshLib
{
namespace ProjectPointOnMesh
{
Element const* getProjectedElement(std::vector<const Element*> const& elements,
                                   Node const& node)
{
    auto is_right_of = [&node](Node const& a, Node const& b) {
        return GeoLib::getOrientationFast(node, a, b) ==
               GeoLib::Orientation::CW;
    };

    for (auto const* e : elements)
    {
        auto const* nodes = e->getNodes();
        if (e->getGeomType() == MeshElemType::TRIANGLE)
        {
            auto const& a = *nodes[0];
            auto const& b = *nodes[1];
            auto const& c = *nodes[2];
            if (!is_right_of(a, b) && !is_right_of(b, c) && !is_right_of(c, a))
            {
                return e;
            }
        }
        else if (e->getGeomType() == MeshElemType::QUAD)
        {
            auto const& a = *nodes[0];
            auto const& b = *nodes[1];
            auto const& c = *nodes[2];
            auto const& d = *nodes[3];
            if (!is_right_of(a, b) && !is_right_of(b, c) &&
                !is_right_of(c, d) && !is_right_of(d, a))
            {
                return e;
            }
        }
    }
    return nullptr;
}

double getElevation(Element const& element, Node const& node)
{
    Eigen::Vector3d const v =
        node.asEigenVector3d() - element.getNode(0)->asEigenVector3d();
    auto const n = FaceRule::getSurfaceNormal(element).normalized();
    return node[2] - n.dot(v) * n[2];
}

}  // namespace ProjectPointOnMesh

}  // end namespace MeshLib
